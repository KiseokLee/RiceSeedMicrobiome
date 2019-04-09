####### Lorenz Curve #########
# Import library
library(ineq)

# Get dataframe of domesticated rice microbiome (Hellinger transformed)
phy.dom <- subset_samples(phy.clean.hell, Type == "Domesticated")
df.dom <- psmelt(phy.dom) %>% group_by(OTU) %>% summarize(total=sum(Abundance))
df.dom <- df.dom$total
df.dom <- df.dom[!df.dom ==0] 

# Plot Lorenz Curve
par(font.axis = 2, font.lab=2)
plot(Lc(df.dom), col='red',lwd=2, xlab='Cumulative percentage of OTUs from the lowest to highest abundance %', ylab='Cumulative share from total abundance %', main='Cumulative share of abundance by OTUs (Lorenz curve)')
ineq(df.dom,type="Gini") 

## Get dataframe of wild rice microbiome (Hellinger transformed)
phy.wil <- subset_samples(phy.clean.hell, Type == "Wild")
df.wil <- psmelt(phy.wil) %>% group_by(OTU) %>% summarize(total=sum(Abundance))
df.wil <- df.wil$total
df.wil <- df.wil[!df.wil ==0]

par(font.axis = 2, font.lab=2)
plot(Lc(df.wil), col='green',lwd=2, xlab='Cumulative percentage of OTUs from the lowest to highest abundance %', ylab='Cumulative share from total abundance %', main='Cumulative share of abundance by OTUs (Lorenz curve)')
ineq(df.wil,type="Gini")

# Let's plot 2 separate graphs into 1 
plot(Lc(df.dom), col='green',lwd=3, xlab='Cumulative percentage of OTUs from the lowest to highest abundance %', ylab='Cumulative share from total abundance %', main='Cumulative share of abundance by OTUs (Lorenz curve)')
lines(Lc(df.wil), col='red',lwd=3)
legend('topleft',pch =c(15), legend=c('Wild','Domesticated'),
       col=c('red','green'))