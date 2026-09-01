library(tidyverse)
library(minpack.lm)
library(broom)

write.csv(shannon_gfm,"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/shannon.gfm_differlevel.csv")
plot_data <- shannon_gfm %>%
  select(Sample,
         Shannon,
         Menhinick_Coef,
         Level,
         Ecosystem.y) %>%
  na.omit() %>% filter(Shannon > 0)


plot_data$Level <- factor(plot_data$Level,
                          levels = c("phylum",
                                     "class",
                                     "order",
                                     "family",
                                     "genus",
                                     "species"))

fit_ecological_model <- function(df, level) {
  if (level %in%
      c("phylum", "class", "order")) {
    model <- lm(Menhinick_Coef ~ Shannon,
                data = df)
    
    type = "Rivet"
    
    
  }
  
  
  
  else if (level == "family") {
    model <- nlsLM(
      Menhinick_Coef ~
        (a * Shannon) / (b + Shannon),
      data = df,
      start = list(
        a = max(df$Menhinick_Coef),
        b = median(df$Shannon)
      )
    )
    
    type = "Redundancy"
    
    
  }
  
  
  
  else if (level %in%
           c("genus", "species")) {
    model <- nlsLM(
      Menhinick_Coef ~
        a * Shannon * exp(-b * Shannon),
      data = df,
      start = list(a = 1,
                   b = 0.1)
    )
    
    type = "Unimodal"
    
    
  }
  
  
  list(model = model,
       type = type)
  
}

curve_data <- plot_data %>%
  group_by(Level) %>%
  group_map(~ {
    level <- as.character(.y$Level)
    
    
    fit <- fit_ecological_model(.x,
                                level)
    
    
    xnew <- seq(min(.x$Shannon, na.rm = TRUE),
                max(.x$Shannon, na.rm = TRUE),
                length.out = 200)
    
    
    pred <- predict(fit$model,
                    newdata = data.frame(Shannon = xnew))
    
    
    data.frame(Shannon = xnew,
               Menhinick_Coef = pred,
               Level = level)
    
    
  }) %>%
  bind_rows()

library(broom)
library(dplyr)


get_model_summary <- function(df, level) {
  fit <- fit_ecological_model(df,
                              level)
  
  
  model <- fit$model
  
  
  # =====================
  # Linear model
  # =====================
  
  if (level %in% c("phylum", "class", "order")) {
    coef_df <- summary(model)$coefficients
    
    r2 <- summary(model)$r.squared
    
    p <- coef_df["Shannon", "Pr(>|t|)"]
    
    
    equation <- sprintf("y = %.3f + %.3f*x",
                        coef(model)[1],
                        coef(model)[2])
    
    
    return(
      tibble(
        Mechanism=fit$type,
        Model="Linear",
        Equation=equation,
        R2=r2,
        P_value=p,
        AIC=AIC(model),
        N=nrow(df)
      )
    )
  }
  
  
  
  # =====================
  # Nonlinear models
  # =====================
  
  
  coef_model <- coef(model)
  
  
  rss <- sum(residuals(model) ^ 2)
  
  
  tss <- sum((df$Menhinick_Coef -
                mean(df$Menhinick_Coef)) ^ 2)
  
  
  r2 <- 1 - rss / tss
  
  
  coef_p <- tryCatch(
    summary(model)$coefficients[, 4],
    
    error = function(e) {
      NA
    }
    
  )
  
  
  
  if (level == "family") {
    equation <- sprintf("y = %.3f*x/(%.3f+x)",
                        coef_model["a"],
                        coef_model["b"])
    
  }
  
  
  if (level %in% c("genus", "species")) {
    equation <- sprintf("y = %.3f*x*exp(-%.3f*x)",
                        coef_model["a"],
                        coef_model["b"])
    
  }
  
  
  
  tibble(
    Mechanism=fit$type,
    Model=class(model)[1],
    Equation=equation,
    R2=r2,
    P_value=min(coef_p,na.rm=TRUE),
    AIC=AIC(model),
    N=nrow(df)
  )
  
}

model_summary <- plot_data %>%
  
  group_by(Level) %>%
  
  group_modify(~{
    
    
    level <- as.character(.y$Level)
    
    
    get_model_summary(
      .x,
      level
    )
    
  }) %>%
  
  ungroup()


model_summary

level_colors <- c(
  "phylum" = "#D73027",
  "class" = "#FC8D59",
  "order" = "#91BFDB",
  "family" = "#4575B4",
  "genus" = "#66BD63",
  "species" = "#762A83"
)
p <- ggplot(plot_data,
            aes(x = Shannon,
                y = Menhinick_Coef,
                color = Level)) +
  
  
  # 原始点
  
  #geom_point(alpha = 0,
   #          size = 2) +
  
  
  # 拟合曲线
  
  geom_line(
    data = curve_data,
    aes(x = Shannon,
        y = Menhinick_Coef,
        color = Level),
    linewidth = 1.3
  ) +
  
  
  scale_color_manual(values = level_colors) +
  
  
  scale_linetype_manual(values = c(
    "Rivet" = "solid",
    "Redundancy" = "dashed",
    "Unimodal" = "dotdash"
  )) +
  
  
  theme_bw() +#facet_wrap(.~Level,scales = "free_x")+
  
  
  theme(
    aspect.ratio = 1,
    
    panel.grid = element_blank(),
    
    legend.position = "right",
    
    axis.text =
      element_text(size = 12,
                   color = "black"),
    
    axis.title =
      element_text(size = 14,
                   face = "bold"),
    
    legend.title =
      element_blank()
    
  ) +
  
  
  labs(x = "Shannon index",
       
       y = "Guild Functional Menhinick")


p


library(tidyverse)
library(minpack.lm)
library(broom)

# ==========================================
# 1. 数据准备
# ==========================================

plot_data <- shannon_gfm %>%
  select(
    Sample,
    Shannon,
    Menhinick_Coef,
    Level,
    Ecosystem.y
  ) %>%
  na.omit() %>%
  filter(Shannon > 0)

plot_data$Level <- factor(
  plot_data$Level,
  levels = c(
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
)

# ==========================================
# 2. 模型拟合函数
# phylum/class/order/family = Redundancy
# genus/species = Unimodal
# ==========================================

fit_ecological_model <- function(df, level) {
  
  # ----------------------------
  # Redundancy hypothesis model
  # y = a*x / (b + x)
  # ----------------------------
  if (level %in% c("phylum", "class", "order", "family")) {
    
    model <- nlsLM(
      Menhinick_Coef ~ (a * Shannon) / (b + Shannon),
      data = df,
      start = list(
        a = max(df$Menhinick_Coef, na.rm = TRUE),
        b = median(df$Shannon, na.rm = TRUE)
      ),
      lower = c(
        a = 0,
        b = 0.0001
      ),
      control = nls.lm.control(maxiter = 1000)
    )
    
    type <- "Redundancy"
  }
  
  # ----------------------------
  # Unimodal pattern
  # y = a*x*exp(-b*x)
  # ----------------------------
  else if (level %in% c("genus", "species")) {
    
    model <- nlsLM(
      Menhinick_Coef ~ a * Shannon * exp(-b * Shannon),
      data = df,
      start = list(
        a = 1,
        b = 0.1
      ),
      lower = c(
        a = 0,
        b = 0.0001
      ),
      control = nls.lm.control(maxiter = 1000)
    )
    
    type <- "Unimodal"
  }
  
  list(
    model = model,
    type = type
  )
}

# ==========================================
# 3. 生成拟合曲线
# ==========================================

curve_data <- plot_data %>%
  group_by(Level) %>%
  group_map(~ {
    
    level <- as.character(.y$Level)
    
    fit <- fit_ecological_model(
      .x,
      level
    )
    
    xnew <- seq(
      min(.x$Shannon, na.rm = TRUE),
      max(.x$Shannon, na.rm = TRUE),
      length.out = 200
    )
    
    pred <- predict(
      fit$model,
      newdata = data.frame(
        Shannon = xnew
      )
    )
    
    data.frame(
      Shannon = xnew,
      Menhinick_Coef = pred,
      Level = level,
      Mechanism = fit$type
    )
    
  }) %>%
  bind_rows()

# ==========================================
# 4. 模型统计结果
# ==========================================

get_model_summary <- function(df, level) {
  
  fit <- fit_ecological_model(
    df,
    level
  )
  
  model <- fit$model
  coef_model <- coef(model)
  
  rss <- sum(residuals(model)^2)
  tss <- sum(
    (
      df$Menhinick_Coef -
        mean(df$Menhinick_Coef)
    )^2
  )
  
  r2 <- 1 - rss / tss
  
  coef_p <- tryCatch(
    summary(model)$coefficients[, 4],
    error = function(e) {
      NA
    }
  )
  
  # ----------------------------
  # Redundancy model equation
  # ----------------------------
  if (level %in% c("phylum", "class", "order", "family")) {
    
    equation <- sprintf(
      "y = %.3f*x/(%.3f+x)",
      coef_model["a"],
      coef_model["b"]
    )
    
    plateau <- coef_model["a"]
    half_saturation <- coef_model["b"]
    
    return(
      tibble(
        Mechanism = fit$type,
        Model = class(model)[1],
        Equation = equation,
        Plateau = plateau,
        Half_saturation = half_saturation,
        Peak_Shannon = NA_real_,
        R2 = r2,
        P_value = min(coef_p, na.rm = TRUE),
        AIC = AIC(model),
        N = nrow(df)
      )
    )
  }
  
  # ----------------------------
  # Unimodal model equation
  # ----------------------------
  if (level %in% c("genus", "species")) {
    
    equation <- sprintf(
      "y = %.3f*x*exp(-%.3f*x)",
      coef_model["a"],
      coef_model["b"]
    )
    
    peak_shannon <- 1 / coef_model["b"]
    peak_y <- (coef_model["a"] / coef_model["b"]) * exp(-1)
    
    return(
      tibble(
        Mechanism = fit$type,
        Model = class(model)[1],
        Equation = equation,
        Plateau = NA_real_,
        Half_saturation = NA_real_,
        Peak_Shannon = peak_shannon,
        Peak_Y = peak_y,
        R2 = r2,
        P_value = min(coef_p, na.rm = TRUE),
        AIC = AIC(model),
        N = nrow(df)
      )
    )
  }
}

model_summary <- plot_data %>%
  group_by(Level) %>%
  group_modify(~ {
    
    level <- as.character(.y$Level)
    
    get_model_summary(
      .x,
      level
    )
    
  }) %>%
  ungroup()

model_summary

# ==========================================
# 5. 颜色
# ==========================================

level_colors <- c(
  "phylum" = "#D73027",
  "class" = "#FC8D59",
  "order" = "#91BFDB",
  "family" = "#4575B4",
  "genus" = "#66BD63",
  "species" = "#762A83"
)

# ==========================================
# 6. 绘图
# ==========================================

 ggplot(
  plot_data,
  aes(
    x = Shannon,
    y = Menhinick_Coef,
    color = Level
  )
) +
  
 # geom_point(
  #  alpha = 0.1,
   # size = 2
  #) +
  
  geom_line(
    data = curve_data,
    aes(
      x = Shannon,
      y = Menhinick_Coef,
      color = Level,
      linetype = Mechanism
    ),
    linewidth = 1.3
  ) +
  
  scale_color_manual(
    values = level_colors
  ) +
  
  scale_linetype_manual(
    values = c(
      "Redundancy" = "solid",
      "Unimodal" = "dotdash"
    )
  ) +
  
  #facet_wrap(
   # . ~ Level,
  #  scales = "free_x"
  #) +
  
  theme_bw() +
  
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.title = element_blank()
  ) +
  
  
  labs(x = "Shannon index",
       
       y = "Guild Functional Menhinick")



