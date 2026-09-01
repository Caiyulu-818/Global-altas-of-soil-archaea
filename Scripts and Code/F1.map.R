# 设置工作目录
setwd("/Users/caiyulu/Desktop/MAGcode/MAG/figure")

# 读取地图数据
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")

# 转换为 Mercator 投影
dat_rob1 <- dat_rob %>%
  # 转换为sf对象(空间数据格式)，指定经纬度列和坐标系(WGS84) 
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  # 转换到罗宾逊投影坐标系
  sf::st_transform(
    crs = sf::st_crs(
      "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
    )
  )

setwd("/Users/caiyulu/Desktop/MAGcode/MAG/figure")
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")
wmap_robin <- spTransform(wmap,CRS("+proj=wintri +datum=WGS84 +no_defs")) #CRS("+proj=robin"))
wmap_df_robin <- fortify(wmap_robin)
# convert to dataframe
wmap_df <- fortify(wmap)

# create a blank ggplot theme
theme_opts <- list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "#97CBFF"),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 22, hjust = .5)
  )
)

arc_sample <-
  read_xlsx(
    "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F1/S1.xlsx"
  )

arc_sample$Lat <- as.numeric(arc_sample$Lat)
arc_sample$Lon <- as.numeric(arc_sample$Lon)


grat <-
  readOGR("ne_110m_graticules_all", layer = "ne_110m_graticules_15")
grat_df <- fortify(grat)

bbox <-
  readOGR("ne_110m_graticules_all", layer = "ne_110m_wgs84_bounding_box")
bbox_df <- fortify(bbox)

grat_robin <-
  spTransform(grat, CRS("+proj=wintri +datum=WGS84 +no_defs"))# CRS("+proj=robin"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <-
  spTransform(bbox,  CRS("+proj=wintri +datum=WGS84 +no_defs"))#CRS("+proj=robin"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

#bbox_mercator <- spTransform(bbox, CRS("+proj=merc"))
#grat_mercator <- spTransform(grat, CRS("+proj=merc"))


s1_sm1 <- arc_sample %>% filter(!Ecosystem == "Sea")
s1_sm1$Ecosystem <- gsub("Glacier", "Wetland", s1_sm1$Ecosystem)
sel <- s1_sm1$Lon > -120 &
  s1_sm1$Lon < -100 &
  s1_sm1$Lat > 0 &
  s1_sm1$Lat < 30

s1_sm1 <- s1_sm1[!sel,]
s1_sm1$Lat <- as.numeric(s1_sm1$Lat)
s1_sm1$Lon <- as.numeric(s1_sm1$Lon)
places_robin_df <- project(cbind(s1_sm1$Lon[], s1_sm1$Lat[]),
                           #proj="+init=ESRI:54030")
                           proj = "+proj=wintri +datum=WGS84 +no_defs")#"+proj=robin")


places_robin_df <- as.data.frame(places_robin_df)
names(places_robin_df) <- c("LONGITUDE", "LATITUDE")

places_robin_df <- places_robin_df %>%
  mutate(
    ecosystem = s1_sm1$Ecosystem,
    sample = s1_sm1$Sample,
    Archaeal_MAGs = s1_sm1$`Archaeal MAGs`,
    Number_of_MAGs = s1_sm1$Genome_num
  ) %>%
  arrange(sample)

write.csv(
  s1_sm1,
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F1/s1_sm1_3671.csv"
)

sample.class <-
  sort(table(places_robin_df$ecosystem), decreasing = TRUE)
leg.lab <- paste0(names(sample.class), " (", sample.class, ")")
yes <- s1_sm1 %>% filter(`Archaeal MAGs` == "Yes")
Archaeal.sample <- sort(table(yes$Ecosystem), decreasing = TRUE)
leg.lab1 <- paste0(names(Archaeal.sample), " (", Archaeal.sample, ")")

# Draw Plot
# 获取世界地图数据
# https://rdrr.io/cran/maps/man/world.html
world <- map_data("world")
p1 <-
  ggplot() + geom_polygon(
    data =world,
    aes(x = long, y = lat, group = group),
    color = "grey20",
    fill = "grey90",
    linewidth = 0.1
  ) +
  theme(panel.background = element_rect(fill = "#eaeef1", colour = "white")) +
  ylim(c(-50, 90)) +
  xlim(c(-180,200))+
  coord_quickmap()
p1
p2 <- p1 +
  geom_point(
    data = na.omit(s1_sm1),
    aes(
      x = Lon,
      y = Lat,
      color = ifelse(`Archaeal MAGs` == "Yes", Ecosystem, "grey"),
      size = ifelse(`Archaeal MAGs` == "Yes", Genome_num, 1),
      shape = ifelse(`Archaeal MAGs` == "Yes", "No", "Yes")
    ),
    # Use 'Yes' and 'No' for shape
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = mycol35,
    name = "Ecosystem\n(Number of samples) of Archaeal MAGs",
    breaks = names(Archaeal.sample),
    labels = leg.lab1
  ) +
  #scale_size_continuous(range = c(0, 490000), name = "Ecosystem\n(Number of MAGs)",
  #   breaks = c(1, 10, 50, 100, 135),
  # limits = c(1, 489732)) +
  scale_shape_manual(values = c("Yes" = 17, "No" = 21), guide = "none") +  # Discrete values for shape
  #theme_void() +
  guides(
    colour = guide_legend(nrow = 3),
    size = guide_legend(title = "Number of Archaeal MAGs", nrow = 1),
    alpha = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.key.size = unit(1, "mm"),
    legend.spacing = unit(1, "mm"),
   # panel.grid =element_blank(),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank()
  )
p2
ggsave(p2,
       filename = "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F1/map.new.pdf",
       height = 8,
       width = 6)


#改变投影方式####

library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)

#分布更换投影方法####
library(sf)
library(rnaturalearth)
library(dplyr)

# 1) 世界多边形（sf）
world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# 2) 先定义投影（示例：Winkel–Tripel）
# 3) 选一个真的不一样的投影（改这一行即可）
# proj <- "+proj=robin"      # Robinson
# proj <- "+proj=wintri"     # Winkel Tripel
proj <- "+proj=eqearth"    # Equal Earth（等积）
#proj <- "+proj=moll"         # Mollweide（等积，椭圆外框，非常好区分）
#proj <- "+proj=wintri +datum=WGS84 +no_defs"
#proj <- "+proj=wintri"
# 3) 在经纬度坐标下做等度网格（20° × 20°），再转投影
bbox_ll  <- st_as_sfc(st_bbox(c(xmin = -180, ymin = -90, xmax = 180, ymax = 90), crs = 4326))
grat_poly_ll <- st_make_grid(bbox_ll, cellsize = c(20, 20), what = "polygons", square = TRUE) |> 
  st_as_sf()

# 4) 全部转到目标投影
world_prj    <- st_transform(world_sf, crs = proj)
grat_poly_prj <- st_transform(grat_poly_ll, crs = proj)
# 2) 你的点数据 -> sf，经纬度列替换成你代码里的列名
# 这里假设 c4 里有 jitter_lon / jitter_lat、Ecosystem、p.abundance 列
# ---- 把 c4 转成 sf，并投影 ----
# 5) 画图：网格是 polygon，可以用 fill 上色
library(ggplot2)

# 找出非数值（例如文本、符号）
s1_sm1_clean<-s1_sm1 %>%
  filter(!grepl("^-?[0-9.]+$", as.character(Lon)) |
           !grepl("^-?[0-9.]+$", as.character(Lat))) %>%
  select(Lon, Lat)

# 找出范围之外的
s1_sm1_clean %>%
  filter(Lon < -180 | Lon > 180 | Lat < -90 | Lat > 90) %>%
  select(Lon, Lat)

# 2) 经纬度 -> sf（注意经度在前、纬度在后）
s1_sm1_sf  <- sf::st_as_sf(s1_sm1_clean, coords = c("Lon", "Lat"), crs = 4326)

# 3) 投影到和底图一致的 Winkel–Tripel（或你定义的 proj）
s1_sm1_prj <- sf::st_transform(s1_sm1_sf, crs = proj)

ggplot() +
  geom_sf(data = grat_poly_prj, fill = "grey89",
          alpha = 0.5, color = NA) + 
  geom_sf(data = world_prj, fill = "grey80", color = "grey99", linewidth = 0.3) +
  geom_sf(data = grat_prj,# fill= "#814662",
          color = "grey70", linewidth = 0.1)+
  geom_sf(
    data = s1_sm1_prj,
    aes(
      color = ifelse(`Archaeal MAGs` == "Yes", Ecosystem, "grey9"),
      size = ifelse(`Archaeal MAGs` == "Yes", Genome_num, 1),
      shape = ifelse(`Archaeal MAGs` == "Yes", "No", "Yes")
    ),
    # Use 'Yes' and 'No' for shape
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = eco_color_map,
    name = "Ecosystem\n(Number of samples) of Archaeal MAGs",
    breaks = names(Archaeal.sample),
    labels = leg.lab1
  ) +
  #scale_size_continuous(range = c(0, 490000), name = "Ecosystem\n(Number of MAGs)",
  #   breaks = c(1, 10, 50, 100, 135),
  # limits = c(1, 489732)) +
  scale_shape_manual(values = c("Yes" = 17, "No" = 21), guide = "none") +  # Discrete values for shape
  #theme_void() +
  guides(
    colour = guide_legend(nrow = 3),
    size = guide_legend(title = "Number of Archaeal MAGs", nrow = 1),
    alpha = "none"
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.box = "vertical",
    legend.box.just = "left",
    legend.key.size = unit(1, "mm"),
    legend.spacing = unit(1, "mm"),
    # panel.grid =element_blank(),
    axis.text =element_blank(),
    axis.ticks =element_blank(),
    axis.title =element_blank()
  )
ggsave(filename = "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F1/map.new_eqearth.pdf",
       height = 10,
       width = 6)

