library(usmap)


usa <- map_data('usa')
ggplot(data=usa, aes(x=long, y=lat, group=group)) + 
  geom_polygon(fill='grey20') + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  coord_fixed(1.3) +
  theme_void()
