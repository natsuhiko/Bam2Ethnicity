GMC <-
function(X,y){
        library(mvtnorm)
        mu=list()
        si=list()
        muall=X
        key = names(table(y))
        nkey=length(key)
        for(i in key){
                mu = c(mu, list(colMeans(X[y%in%i,])))
                muall[y%in%i,]=rep(1,sum(y%in%i,na.rm=T))%*%t(colMeans(X[y%in%i,]))
                #si = c(si, list(cov(X[y%in%i,])))
        }
        names(mu)=key
        keypair=NULL
        for(i in 1:(nkey-1)){
                for(j in (i+1):nkey){
                        keypair=c(keypair, paste(key[i],key[j],sep=":"))
                        mu = c(mu, list((mu[[key[i]]]+mu[[key[j]]])/2))
                }
        }
        names(mu)=c(key,keypair)
        for(i in names(mu)){ si = c(si, list(cov(X-muall)*40)) }
        Z=array(0,c(sum(is.na(y)),length(mu)))
        k=1
        priorss=rep(1,length(mu))
        for(i in length(si)){si[[i]]=si[[i]]*priorss[i]}
        for(i in seq(nrow(X))[is.na(y)]){
                for(j in seq(length(mu))){
                        Z[k,j] = dmvnorm(X[i,],mu[[j]],si[[j]],log=T)
                }
                k=k+1
        }
        colnames(Z)=names(mu)
        invisible(exp(Z)/rowSums(exp(Z)))
}


spop=c(scan("/nfs/team205/nk5/Applications/EthnicityEst/spop.txt",""),NA)
evec=read.table(commandArgs()[5])
print(dim(evec))
print(head(evec))
P=GMC(evec[!spop%in%"AMR",2:4],spop[!spop%in%"AMR"])
write.table(P[nrow(P),,drop=F],file=commandArgs()[6],row=F,col=T,sep="\t",quote=F)

