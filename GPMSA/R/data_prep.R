## prepare sim data
## eta_n x m matrix where m is the number of unique input settings
## and eta_n is the length of each output

PDF=F

simout_path <- '../../ExPR1/out/ott_51.05/'
files <- list.files(simout_path, pattern = "csv")
nfile <- length(files)

m <- nfile
nlag <- 1
sim_native <- matrix(NA, nrow = 90*nlag, ncol = m)
# x_native <- matrix(NA, nrow = m, ncol = 2)
x_native <- matrix(NA, nrow = m, ncol = 1)
  
for (i in 1:nfile){
  # read data
  d <- read.csv(paste0(simout_path, files[i]), stringsAsFactors = F)
  sim_native[,i] <- d$euvac
  # extract input from file name
  # x_native[i,1] <- as.numeric(unlist(strsplit(files[i], '_'))[2])
  x_native[i,1] <- as.numeric(unlist(strsplit(unlist(strsplit(files[i], '_'))[2], '.csv')))
  # x_native[i,2] <- as.numeric(unlist(strsplit(files[i], '_'))[3])
  # x_native[i,3] <- as.numeric(unlist(strsplit(unlist(strsplit(files[i], '_'))[4], '.csv')))
}

if(PDF) jpeg('../plots/design.jpeg', width=880, height=880)
# pairs(x_native, labels = c(expression(eta), expression(alpha)), gap=0, cex.labels = 6)
par(mfrow=c(1,1), mar=c(8,8,3,3))
plot(x_native, xlab = expression(eta), ylab = expression(alpha), cex.lab=3, cex = 5, pch=19, col = 'grey')
abline(v=seq(0,0.8, by = 0.1), col = 'grey')
abline(h=seq(0.8,3, by = 0.5), col = 'grey')
if(PDF) dev.off()


## read observed data

obs <- read.csv('../../Data/csv/20150311_ott_abs.csv', stringsAsFactors = F)
# obs <- ifelse(obs[,2]<0, 0, obs[,2])
obs <- obs[,2]

## a simple plot of obs and sim

if(PDF) jpeg('../plots/sim-obs.jpeg', width = 1000, height = 600)
par(mfrow=c(1,1))
plot(NA, xlim = c(0,90), ylim = c(0,5), axes = F, ylab = "absorption (in dB)", xlab = 'time (seconds)')
matplot(sim_native, col = 'grey', type = "l", add = T)
points(obs, col = 'red', lty = 2, lwd = 2, pch = 19)
axis(1)
axis(2)
box()
legend(60,4,c('simulated','Actual'), col=c('grey','red'),lty = c(3,NA), cex = 2, pch = c(NA, 19))
if(PDF) dev.off()

write.table(x_native, file="design_ott.txt", quote = F, row.names = F, col.names = F)
write.table(matrix(obs, ncol = 1), file="obs_ott.txt", quote = F, row.names = F, col.names = F)
write.table(sim_native, file="sim_ott.txt", quote = F, row.names = F, col.names = F)
