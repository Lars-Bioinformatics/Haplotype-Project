library(viridis)
library(ggplot2)
library(dplyr)

map.world <- map_data('world')
map.world <- map.world %>% filter(region != "Antarctica")
map.world$region <- toupper(map.world$region)
#map.world <- fortify(map.world)

# European countries
europe <- toupper(c("Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina", 
                    "Bulgaria", "Croatia", "Czech Republic", "Denmark", "Estonia", "Finland", "France", 
                    "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", 
                    "Liechtenstein", "Lithuania", "Luxembourg", "Malta", "Moldova", "Monaco", 
                    "Netherlands", "Norway", "Poland", "Portugal", "Romania", "San Marino", 
                    "Serbia and Montenegro", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", 
                    "Ukraine", "UK"))

russia = "Russia"

northAmerica = toupper(c("USA", "Canada", "Mexico", "Greenland"))

if (sum(unique(brca_data$Country) %in% europe) == length(unique(brca_data$Country))){
    map <- filter(map.world, region %in% europe)
    # remove svalbard from Norway
    map <- filter(map, lat < 73)
} else if (sum(unique(brca_data$Country) %in% c(europe, northAmerica)) == length(unique(brca_data$Country))){
    map <- filter(map.world, long < 50, long > -250, lat > 10, lat < 150)#, region %in% c(europe, northAmerica))
    #map <- filter(map.world, long < 50, long > -250, region %in% c(europe, northAmerica))
} else {
    map <- map.world
}

map %>% group_by(region) %>% summarise() %>% print(n = Inf)

# Keep founder groups
brca_data2 = subset(brca_data, cluster_groups %in% filter(count(brca_data,cluster_groups), n>4)$cluster_groups)

hf <- suppressWarnings(left_join(map, brca_data2, by = c('region' = 'Country')) %>% filter(Mut1HGVS == mut))
#hf$Mut1HGVS = as.character(hf$Mut1HGVS)

# countries <- map.world %>% filter(region %in% brca_data$Country)
# country_centriods <- do.call(rbind, sapply(unique(countries$region), function(country) maps:::apply.polygon(subset(countries, region == country), maps:::centroid.polygon)))
# country_centriods = setNames(as.data.frame(country_centriods), c("x","y"))

countries <- read.table("/Users/lars/Desktop/country_centroids.txt", header = T, sep = "\t")
countries$name <- toupper(countries$name)

country_centriods = filter(countries, name %in% unique(brca_data$Country))
country_centriods

# arrows_straight = data.frame(xstart=country_centriods$longitude[seq(1,5,2)], ystart=country_centriods$latitude[seq(1,5,2)], xend=country_centriods$longitude[seq(2,6,2)], yend=country_centriods$latitude[seq(2,6,2)])
# arrows_curve = data.frame(xstart=country_centriods$longitude[seq(2,5,2)], ystart=country_centriods$latitude[seq(2,5,2)], xend=country_centriods$longitude[seq(3,6,2)], yend=country_centriods$latitude[seq(3,6,2)])

for (i in unique(brca_data2$cluster_groups)){
    brca_data.temp = subset(brca_data, cluster_groups == i)
    
    if (sum(unique(brca_data.temp$Country) %in% europe) == length(unique(brca_data.temp$Country))){
        map <- filter(map.world, region %in% europe)
        # remove svalbard from Norway
        map <- filter(map, lat < 73)
        # Heigh and with of plot
        width = 13; height = 8; scale = 0.6
    } else if (sum(unique(brca_data.temp$Country) %in% c(europe, northAmerica)) == length(unique(brca_data.temp$Country))){
        map <- filter(map.world, long < 50, long > -250, lat > 10, lat < 150)#, region %in% c(europe, northAmerica))
        #map <- filter(map.world, long < 50, long > -250, region %in% c(europe, northAmerica))
        width = 14; height = 5; scale = 0.9
    } else {
        map <- map.world
    }
    
    hf <- suppressWarnings(left_join(map, brca_data.temp, by = c('region' = 'Country')) %>% filter(Mut1HGVS == mut))
    country_centriods = filter(countries, name %in% unique(brca_data.temp$Country))
    
    if (length(unique(brca_data.temp)) == 1){ # DENMARK
        arrows_straight = data.frame()
        arrows_curve = data.frame()
    } else { # SPAIN
        arrows_straight = data.frame()
        arrows_curve = data.frame(xstart=country_centriods$longitude[2], xend=country_centriods$longitude[1], ystart=country_centriods$latitude[2], yend=country_centriods$latitude[1])
    }
    
    cols <- colors_plot(brca_data)
    cols <- cols[names(cols) %in% brca_data.temp$Country]
    p <- ggplot() + #ggplot(data = hf, aes(x=long, y=lat, group=group)) +
        geom_map(data = map, map = map, aes(map_id=region), fill="lightgrey", color="grey") +
        geom_polygon(data = hf, aes(x=long, y=lat, group=group, fill=region)) +
        #geom_point(data = country_centriods, aes(x=longitude, y=latitude), color = "darkblue", size = 1) +
        #geom_line(data = country_centriods, aes(x=longitude, y=latitude), arrow=arrow(), lineend = "round", linejoin = "round") +
        coord_fixed(1.3) +
        # coord_equal(1.3) +
        expand_limits(x = map$long, y = map$lat) +
        scale_fill_manual(values = cols) + 
        # theme_light()
        # theme_minimal()
        theme_bw()
    if (nrow(arrows_curve)>0) {
        p <- p + geom_curve(data = arrows_curve, aes(x=xstart, y=ystart, xend=xend, yend=yend), arrow=arrow(length = unit(0.3, "cm")), lineend = "round")#, linejoin = "round") +
    }
    if (nrow(arrows_straight)>0) {
        p <- p + geom_segment(data = arrows_straight, aes(x=xstart, y=ystart, xend=xend, yend=yend), arrow=arrow(length = unit(0.3, "cm")), lineend = "round", linejoin = "round")
    }
    p
    # ggsave(filename = paste0("cache/BRCA2/migration/c.7617+1G_A/BRCA2-c.7617+1G_A-Migration-group", i, ".pdf"), plot = p, width = width, height = height, scale = scale)
    ggsave(filename = paste0("cache/BRCA2/migration/c.7617+1G_A/BRCA2-c.7617+1G_A-Migration-group", i, ".png"), plot = p, width = width, height = height, scale = scale)
}
