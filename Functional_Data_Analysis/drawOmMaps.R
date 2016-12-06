library(rworldmap)
newmap <- getMap(resolution = "low")
plot(newmap,xlim = c(-150, -50), ylim = c(30, 80), asp = 1)
plot(newmap)
points(x=fbusi$longitude,y=fbusi$latitude, pch=19, col='red' , cex=0.01)
points(-10,0, lwd=20, col='red')

library(ggmap)
map=get_map('Canada', zoom = 3)
ggmap(map) + geom_point(data=w, aes(x=-W.longitude, y=N.latitude, col=region), size=3)

