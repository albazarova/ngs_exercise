#Summary: I develop a Bayesian methodology to identify patients whose malignant isoform counts increased independently of their non-malignant isoform counts
#I build a model mu_i+eps for counts, where mu_i is the true number of counts for the i-th patient
#eps is noise with zero mean and infinite variance, assumed to follow Student distribution with 2 degrees of freedom, 0 location and greater than 1 scale parameter
#I heuristically find the scale and degrees of freedom parameters by observing the difference between measurements of the counts and sample mean of these measurements for each patient
#I discover that noise follows different distribution for patients 1 to 2160 and the remaining ones
#I identify the noise parameter for patients 1:2160 non-malignant isoform (i), 1:2160 malignant isoform (ii), remaining patients non-malignant isoform (iii), remaining patients malignant isoform (iv), thereby divideing patients into two groups 

#Further in a Bayesian MCMC algorithm I attempt to effectively subtract the noise from the measurement in order to obtain the true count values for each of the patients along with its confidences and distributions over these confidences
#The acceptance rate turned out quite low which indicates that the Metropolis-Hastings algorithm within the MCMC sampler should be improved, firstly by choosing a more appropriate proposal
#That was also partly due to the number of missing values in the second group of patients
#Partly because of that there were much less patients with significant changes identified in the second group

#I dealt with the missing values by not including them into the likelihood
#Negative values were excluded, as well as abnormally large ones, displayed as Inf by the software

#Patients with significant changes in malignant isoform were identified by investigating the resulting distributions of count parameters obtained by MCMC algorithm
#I computed the difference in mean values of counts for malignant and non-malignant isoform counts in December and in August, as well as the probability of these counts coming from the same distribution
#Based on the latter criteria I choose the ones with significant changes and filter out those ones who had a change in malignant isoform but not in non-malignant

#The obtained solutions (i.e. patients with significant changes in malignant isoform) were written into the files part1.csv and part2.csv (please see the repository)

#Overall the developed methodology is certainly capable of dealing with the given problem, but would require further work in order to produce robust results


library(matrixStats)

#load data
data<-read.csv("https://raw.githubusercontent.com/albazarova/ngs_exercise/master/exercise.csv")

#identify patients whi have non-missing August measurement of malignancy (otherwise there is nothing to compare with)

a<-array(0,dim(data)[1])

ind=c();

for (i in 1:dim(data)[1]){
  if (sum(is.na(data[i,1:4]))<4){
    ind=c(ind,i)
  }
}

#to make conclusions about the sequencing noise makes sense to subtract mean observation from each of the patients

noise_ma<-data[ind,1:4]-rowMeans(data[ind,1:4],na.rm=TRUE)
noise_md<-data[ind,5:8]-rowMeans(data[ind,5:8],na.rm=TRUE)
noise_ta<-data[ind,9:12]-rowMeans(data[ind,9:12],na.rm=TRUE)
noise_td<-data[ind,13:16]-rowMeans(data[ind,13:16],na.rm=TRUE)

#noise is obviously additive with mean 0, also observe that after patient 2160 its variance becomes significantly smaller and therefore makes sense to model patients 1:2160 and 2160:14431 separately
#the above is obvious from the plot, e.g.

plot(unlist(noise_ma[,1]),type='l')
plot(unlist(noise_ta[,1]),type='l')

#also it is obvious that for malignant and non-malignant isoforms noise is different and therefore have also to be modeled separately

#further taking into account the heavy tails of the distribution we model noise as a generalised Student distribution with 2 degrees of freedom, 0 location parameter and variable scale parameter, i.e.
#scale=130 for first 2160 patients when modelling malignant isoform
#scale=500 for first 2160 patients when modelling non-malignant isoform
#scale=30 for the remaining patients when modelling malignant isoform
#scale=800 for the remaining patients when modelling non-malignant isoform
#all these functions seem to be a good fit for noise, e.g.

par(mfrow=c(1,2))
 hist(noise_ma[1:2160,1],100)
 hist(rt(10000,2)*130,100)
 
#The following function computes Student density with given parameters at a given point
 
Stu<-function(x, mu,sig,nu ){

 r<-gamma((nu+1)/2)*((1+(1/(nu*sig))*(x-mu)^2)^(-(nu+1)/2))/(gamma(nu/2)*sqrt(pi*nu*sig));
 
 r
}

#We model each measurement as mu+eps, where mu is the true count value, same for all 4 observation, and eps is the Student noise with appropriate parameters (please see above)
#Our aim is to effectively subtract noise from the measurements in order to obtain the true count numbers and their confidences for each of the parients
#We achieve that by utilising Bayesian inference with Metropolice-Hastings sampler
#Ideally the parameters of the Student distribution should be inferred, but here to save time we assume they are known. The algorithm however is written in such a way that it will be easy to generalise it to the case when noise parameters are unknown
#Also convergence rate can be improved by setting priors on mu parameters

#The following code is for implementation of the above methodology for August non-malignant isoforms

N<-2159; #number of patients
it<-10000; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(500^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,12)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=3000; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){

  print(i) #to know which iteration is running

  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate

  ratma=0;
ratm=0;
}


for (j in 1:N){ #for each patient a separate proposal

ratm=ratm+1;

mup[j]=rnorm(1,mu[j,i],s); #sample proposal

while(mup[j]<0){
  mup[j]=rnorm(1,mu[j,i],s);
}

if (j==1){#computing likelihood
  li1=sum(log(Stu(data[1:N,9],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,10],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,11],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,12],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
lip=sum(log(Stu(data[1:N,9],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,10],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,11],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,12],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
} else{
  li1=sum(log(Stu(data[1:N,9],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,10],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,11],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,12],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
lip=sum(log(Stu(data[1:N,9],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,10],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,11],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,12],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
}
alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals

if (alpha>runif(1)){
  mu[j,i+1]=mup[j];
ratma=ratma+1;
} else {
  mu[j,i+1]=mu[j,i];
}

}



}

#mu is the result of the MCMC simulation; ideally the MCMC chain has to be tested for convergence via GelmanRubin statistics, but here we skip it
#we write the obtained chain into a file
#write.table(mu,"datat_1.txt",sep="\t");

#The following code is for implementation of the above methodology for December non-malignant isoforms for fist 2160

N<-2158 #negative values observed at position 2160 and 128, they have to be eliminated

data<-rbind(data[1:127,],data[129:2159,],data[2161:22089,]);

it<-10000; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(500^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,12)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=3000; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[1:N,13],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,14],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,15],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,16],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,13],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,14],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,15],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,16],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[1:N,13],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,14],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,15],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,16],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,13],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,14],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,15],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,16],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation, we write it into a file
#write.table(mu,"datat_2.txt",sep="\t");

#The following code is for implementation of the above methodology for August malignant isoforms

N<-2160; #number of patients
it<-7000; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(130^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,10)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=500; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[1:N,1],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,2],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,3],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,4],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,1],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,2],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,3],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,4],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[1:N,1],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,2],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,3],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,4],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,1],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,2],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,3],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,4],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation; ideally the MCMC chain has to be tested for convergence via GelmanRubin statistics, but here we skip it
#we write the obtained chain into a file
#write.table(mu,"datam_1.txt",sep="\t");

#The following code is for implementation of the above methodology for December non-malignant isoforms for fist 2160

N<-2160 

it<-7000; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(130^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,12)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=500; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[1:N,5],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,6],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,7],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,8],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,5],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,6],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,7],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,8],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[1:N,5],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,6],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,7],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,8],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,5],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,6],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,7],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,8],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation, we write it into a file
#write.table(mu,"datam_2.txt",sep="\t"); this is a rather large file


#We are now going to perform the analysis of the difference between counts in December and August
#For that we will compute the difference between mean values of the inferred mus (December - August)
#We are going to compute the probability of the count values in December and August coming from the same distribution, by integrating intersecting area of their distribution
#For both malignant and non-malignant isoforms we will filter out those which exhibited difference in mean values and the probability of the counts coming from the same distribution is below some parameter value (e.g. 0.05)
#We will then identify those patients for whom significant difference was detected in malignant isoform, but not in non-malignant

pr<-0.1 #parameter for probability

#s1<-read.table("datat_1.txt"); #August t
#s<-read.table("datat_2.txt"); #December t
N=2158; #because of two instances of negative values
o=10000;
b=o/2;

box=1000;

mea<-array(0,N); #December-August mean
prob<-array(0,N); #Probabilities

#the following code finds the difference in the means and writes it into mea[]
#and the probabilities of coming from the same distribution and writes them into prob[]
for (i in 1:127){
  
  print(i)
mea[i]=mean(as.numeric(s[i,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i,(b+1):o]),na.rm=TRUE);

hist1<-hist(unlist(s1[i,(b+1):o]),box,plot=FALSE); #distribution of the mus in August
hist2<-hist(unlist(s[i,(b+1):o]),box,plot=FALSE); #December
h<-union(hist1$mids,hist2$mids)

if (min(length(hist1$counts),length(hist2$counts))>1){
c10<-approx(hist1$mids,hist1$counts,h)
c20<-approx(hist2$mids,hist2$counts,h)
c10$y[which(is.na(c10$y)==TRUE)]=0;
c20$y[which(is.na(c20$y)==TRUE)]=0;
prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
} else{
 
  if (max(length(hist1$counts),length(hist2$counts))==1) {
    prob[i]=0;
  } else{
    
    if (length(hist1$counts)==1){
      
      c20<-approx(hist2$mids,hist2$counts,h)
      
      c20$y[which(is.na(c20$y)==TRUE)]=0;
      
      prob[i]=c20$y[hist1$mids];
    } else{
      
      c10<-approx(hist1$mids,hist1$counts,h)
      c10$y[which(is.na(c10$y)==TRUE)]=0;
      
      prob[i]=c10$y[hist2$mids];
    }
    
  }
}
}

for (i in 128:N){
  
  print(i)
  
  mea[i]=mean(as.numeric(s[i,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i+1,(b+1):o]),na.rm=TRUE);
  hist1<-hist(unlist(s1[i+1,(b+1):o]),box,plot=FALSE);
  hist2<-hist(unlist(s[i,(b+1):o]),box,plot=FALSE);
  h<-union(hist1$mids,hist2$mids)
  if (min(length(hist1$counts),length(hist2$counts))>1){
    c10<-approx(hist1$mids,hist1$counts,h)
    c20<-approx(hist2$mids,hist2$counts,h)
    c10$y[which(is.na(c10$y)==TRUE)]=0;
    c20$y[which(is.na(c20$y)==TRUE)]=0;
    prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
  } else{
    
    if (max(length(hist1$counts),length(hist2$counts))==1) {
      prob[i]=0;
    } else{
      
      if (length(hist1$counts)==1){
        
        c20<-approx(hist2$mids,hist2$counts,h)
        
        c20$y[which(is.na(c20$y)==TRUE)]=0;
        
        prob[i]=c20$y[hist1$mids];
      } else{
        
        c10<-approx(hist1$mids,hist1$counts,h)
        c10$y[which(is.na(c10$y)==TRUE)]=0;
        
        prob[i]=c10$y[hist2$mids];
      }
      
    }
  }
}

prob_t<-prob

mea_t<-mea
 
 
s1<-read.table("/Users/bazarova/Dropbox/datam_1.txt"); #August t
s<-read.table("/Users/bazarova/Dropbox/datam_2.txt"); #December t
N=2158; #because of two instances of negative values
o=7000;
b=o/2;

box=700;

mea<-array(0,N); #December-August mean
prob<-array(0,N); #Probabilities

#the following code finds the difference in the means and writes it into mea[]
#and the probabilities of coming from the same distribution and writes them into prob[]
for (i in 1:127){
  
  print(i)
  mea[i]=mean(as.numeric(s[i,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i,(b+1):o]),na.rm=TRUE);
  
  hist1<-hist(unlist(s1[i,(b+1):o]),box,plot=FALSE); #distribution of the mus in August
  hist2<-hist(unlist(s[i,(b+1):o]),box,plot=FALSE); #December
  h<-union(hist1$mids,hist2$mids)
  
  if (min(length(hist1$counts),length(hist2$counts))>1){
    c10<-approx(hist1$mids,hist1$counts,h)
    c20<-approx(hist2$mids,hist2$counts,h)
    c10$y[which(is.na(c10$y)==TRUE)]=0;
    c20$y[which(is.na(c20$y)==TRUE)]=0;
    prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
  } else{
    
    if (max(length(hist1$counts),length(hist2$counts))==1) {
      prob[i]=0;
    } else{
      
      if (length(hist1$counts)==1){
        
        c20<-approx(hist2$mids,hist2$counts,h)
        
        c20$y[which(is.na(c20$y)==TRUE)]=0;
        
        prob[i]=c20$y[hist1$mids];
      } else{
        
        c10<-approx(hist1$mids,hist1$counts,h)
        c10$y[which(is.na(c10$y)==TRUE)]=0;
        
        prob[i]=c10$y[hist2$mids];
      }
      
    }
  }
}

for (i in 128:N){
  
  print(i)
  
  mea[i]=mean(as.numeric(s[i+1,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i+1,(b+1):o]),na.rm=TRUE);
  hist1<-hist(unlist(s1[i+1,(b+1):o]),box,plot=FALSE);
  hist2<-hist(unlist(s[i+1,(b+1):o]),box,plot=FALSE);
  h<-union(hist1$mids,hist2$mids)
  if (min(length(hist1$counts),length(hist2$counts))>1){
    c10<-approx(hist1$mids,hist1$counts,h)
    c20<-approx(hist2$mids,hist2$counts,h)
    c10$y[which(is.na(c10$y)==TRUE)]=0;
    c20$y[which(is.na(c20$y)==TRUE)]=0;
    prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
  } else{
    
    if (max(length(hist1$counts),length(hist2$counts))==1) {
      prob[i]=0;
    } else{
      
      if (length(hist1$counts)==1){
        
        c20<-approx(hist2$mids,hist2$counts,h)
        
        c20$y[which(is.na(c20$y)==TRUE)]=0;
        
        prob[i]=c20$y[hist1$mids];
      } else{
        
        c10<-approx(hist1$mids,hist1$counts,h)
        c10$y[which(is.na(c10$y)==TRUE)]=0;
        
        prob[i]=c10$y[hist2$mids];
      }
      
    }
  }
}

prob_m<-prob

mea_m<-mea


#indices of counts which changed from August to December
ind_t=intersect(which((prob_t<pr)&(prob_t>0)),which(abs(mea_t)>0));
ind_m=intersect(which((prob_m<0.1)&(prob_m>0)),which(abs(mea_m)>0));

ind_1=setdiff(ind_m,ind_t) #the patients whose counts in malignant isoform significantly changed

ind_1[ind_1>127]=ind_1[ind_1>127]+1

#write.csv(ind_1,"part1.csv"); write into a file


#The following is applied for the remaining patients who did not have malignancy counts missing in August, approximately 12000 patients
#Everything works in a similar way, however there are a lot of missing data. We skip them by not computing likelihood in those points
#The acceptance rate is very poor suggesting that more work should be done to ensure the proposal is chosen appropriately
#The computational complexity also grows and therefore it would be better to transform it into a C++ code

#analysis for non-malignant counts August

#First, get rid of the Inf values

ind1=intersect(which((data[ind[2161:14431],9]!=Inf)|(is.na(data[ind[2161:14431],9])==TRUE)),which((data[ind[2161:14431],10]!=Inf)|(is.na(data[ind[2161:14431],10])==TRUE)));
for (i in 11:16){
ind1=intersect(ind1,which((data[ind[2161:14431],i]!=Inf)|(is.na(data[ind[2161:14431],i])==TRUE)));
}
ind2=ind[2161:14431];

N=length(ind1);#number of patients


it<-570; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(800^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,10)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=10000; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[ind2[ind1],9],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],10],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],11],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],12],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind2[ind1],9],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],10],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],11],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],12],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[ind2[ind1],9],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],10],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],11],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],12],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[1:N,9],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,10],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,11],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[1:N,12],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation; ideally the MCMC chain has to be tested for convergence via GelmanRubin statistics, but here we skip it
#we write the obtained chain into a file
#write.table(mu,"datat_1_2.txt",sep="\t");

#The following code is for implementation of the above methodology for December non-malignant isoforms for the remaining patients


it<-570; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(800^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,12)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=10000; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[ind2[ind1],13],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],14],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],15],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],16],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind2[ind1],13],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],14],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],15],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],16],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[ind2[ind1],13],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],14],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],15],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],16],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind2[ind1],13],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],14],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],15],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind2[ind1],16],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation, we write it into a file
#write.table(mu,"datat_2_2.txt",sep="\t");

#The following code is for implementation of the above methodology for August malignant isoforms

N<-14431-2160; #number of patients
it<-650; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(30^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,10)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=100; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[ind[2161:2160+N],1],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],2],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],3],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],4],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind[2161:2160+N],1],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],2],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],3],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],4],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[ind[2161:2160+N],1],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],2],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],3],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],4],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind[2161:2160+N],1],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],2],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],3],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],4],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation; ideally the MCMC chain has to be tested for convergence via GelmanRubin statistics, but here we skip it
#we write the obtained chain into a file
#write.table(mu,"datam_1_2.txt",sep="\t");

#The following code is for implementation of the above methodology for December non-malignant isoforms for fist 2160

N<-14431-2160

it<-650; #number of iterations

nu<-2*array(1,it); #degrees of freedom for Student dist
sig<-(30^2)*array(1,it); #scale parameter for Student dist
mu<-matrix(0,N+1,it); #matrux of true count value parameters

mu[,1]<-exp(runif(N,4,12)); #initialising the count values

ratm=0; #for computing the acceptance rate
ratma=0; #for computing the acceptance rate
s=100; #parameter of the random walk proposal, chosen heuristically, ideakky should be tuned during the algorithm

for (i in 1:it){
  
  print(i) #to know which iteration is running
  
  mup<-array(1,N); #initialise the proposal
  
  if (it%%10==0){ #check the acceptance rate every 10 iterations
    
    print(ratma/ratm); #print acceptance rate
    
    ratma=0;
    ratm=0;
  }
  
  
  for (j in 1:N){ #for each patient a separate proposal
    
    ratm=ratm+1;
    
    mup[j]=rnorm(1,mu[j,i],s); #sample proposal
    
    while(mup[j]<0){
      mup[j]=rnorm(1,mu[j,i],s);
    }
    
    if (j==1){#computing likelihood
      li1=sum(log(Stu(data[ind[2161:2160+N],5],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],6],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],7],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],8],mu[1:N,i],sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind[2161:2160+N],5],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],6],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],7],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],8],c(mup[j],mu[2:N,i+1]),sig[i+1],nu[i+1])),na.rm=TRUE);
    } else{
      li1=sum(log(Stu(data[ind[2161:2160+N],5],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],6],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],7],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],8],c(mu[1:(j-1),i+1],mu[j:N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
      lip=sum(log(Stu(data[ind[2161:2160+N],5],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],6],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],7],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE)+sum(log(Stu(data[ind[2161:2160+N],8],c(mu[1:(j-1),i+1],mup[j],mu[(j+1):N,i]),sig[i+1],nu[i+1])),na.rm=TRUE);
    }
    alpha=exp(lip-li1)*((1-dnorm(mu[j,i]/s))/(1-dnorm(mup[j]/s))); #likelihood times the ratio of conditional proposals
    
    if (alpha>runif(1)){
      mu[j,i+1]=mup[j];
      ratma=ratma+1;
    } else {
      mu[j,i+1]=mu[j,i];
    }
    
  }
  
  
  
}

#mu is the result of the MCMC simulation, we write it into a file
#write.table(mu,"datam_2_2.txt",sep="\t"); this is a rather large file


#We are now going to perform the analysis of the difference between counts in December and August
#For that we will compute the difference between mean values of the inferred mus (December - August)
#We are going to compute the probability of the count values in December and August coming from the same distribution, by integrating intersecting area of their distribution
#For both malignant and non-malignant isoforms we will filter out those which exhibited difference in mean values and the probability of the counts coming from the same distribution is below some parameter value (e.g. 0.05)
#We will then identify those patients for whom significant difference was detected in malignant isoform, but not in non-malignant

pr<-0.1 #parameter for probability

s1<-read.table("/Users/bazarova/OneDrive/MATLAB/datam_1_2.txt"); #August t
s<-read.table("/Users/bazarova/Dropbox/datam_2_2.txt"); #December t
N=length(ind1); 
o=570;
b=o/2;

box=50;

mea<-array(0,N); #December-August mean
prob<-array(0,N); #Probabilities

#the following code finds the difference in the means and writes it into mea[]
#and the probabilities of coming from the same distribution and writes them into prob[]
for (i in 1:N){
  
  print(i)
  mea[i]=mean(as.numeric(s[i,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i,(b+1):o]),na.rm=TRUE);
  
  hist1<-hist(unlist(s1[i,(b+1):o]),box,plot=FALSE); #distribution of the mus in August
  hist2<-hist(unlist(s[i,(b+1):o]),box,plot=FALSE); #December
  h<-union(hist1$mids,hist2$mids)
  
  if (min(length(hist1$counts),length(hist2$counts))>1){
    c10<-approx(hist1$mids,hist1$counts,h)
    c20<-approx(hist2$mids,hist2$counts,h)
    c10$y[which(is.na(c10$y)==TRUE)]=0;
    c20$y[which(is.na(c20$y)==TRUE)]=0;
    prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
  } else{
    
    if (max(length(hist1$counts),length(hist2$counts))==1) {
      prob[i]=0;
    } else{
      
      if (length(hist1$counts)==1){
        
        c20<-approx(hist2$mids,hist2$counts,h)
        
        c20$y[which(is.na(c20$y)==TRUE)]=0;
        
        prob[i]=c20$y[hist1$mids];
      } else{
        
        c10<-approx(hist1$mids,hist1$counts,h)
        c10$y[which(is.na(c10$y)==TRUE)]=0;
        
        prob[i]=c10$y[hist2$mids];
      }
      
    }
  }
}


prob_t_2<-prob

mea_t_2<-mea


s1<-read.table("/Users/bazarova/Dropbox/datat_1_2.txt"); #August t
s<-read.table("/Users/bazarova/Dropbox/datat_2_2.txt"); #December t
N=14431-2160; 
o=650;
b=o/2;

box=50;

mea<-array(0,N); #December-August mean
prob<-array(0,N); #Probabilities

#the following code finds the difference in the means and writes it into mea[]
#and the probabilities of coming from the same distribution and writes them into prob[]
for (i in 1:N){
  
  print(i)
  mea[i]=mean(as.numeric(s[i,(b+1):o]),na.rm=TRUE)-mean(as.numeric(s1[i,(b+1):o]),na.rm=TRUE);
  
  hist1<-hist(unlist(s1[i,(b+1):o]),box,plot=FALSE); #distribution of the mus in August
  hist2<-hist(unlist(s[i,(b+1):o]),box,plot=FALSE); #December
  h<-union(hist1$mids,hist2$mids)
  
  if (min(length(hist1$counts),length(hist2$counts))>1){
    c10<-approx(hist1$mids,hist1$counts,h)
    c20<-approx(hist2$mids,hist2$counts,h)
    c10$y[which(is.na(c10$y)==TRUE)]=0;
    c20$y[which(is.na(c20$y)==TRUE)]=0;
    prob[i]=sum(rowMins(cbind(c10$y,c20$y)),na.rm=TRUE)/mean(c(sum(c10$y),sum(c20$y)));
  } else{
    
    if (max(length(hist1$counts),length(hist2$counts))==1) {
      prob[i]=0;
    } else{
      
      if (length(hist1$counts)==1){
        
        c20<-approx(hist2$mids,hist2$counts,h)
        
        c20$y[which(is.na(c20$y)==TRUE)]=0;
        
        prob[i]=c20$y[hist1$mids];
      } else{
        
        c10<-approx(hist1$mids,hist1$counts,h)
        c10$y[which(is.na(c10$y)==TRUE)]=0;
        
        prob[i]=c10$y[hist2$mids];
      }
      
    }
  }
}



prob_m_2<-prob

mea_m_2<-mea


#indices of counts which changed from August to December
ind_t_2=intersect(which((prob_t_2<pr)&(prob_t_2>0)),which(abs(mea_t_2)>0));

ind_t_2=ind2[ind1[ind_t_2]]

ind_m_2=intersect(which((prob_m_2<0.1)&(prob_m_2>0)),which(abs(mea_m_2)>0));

ind_m_2=ind[2160+ind_m_2]

ind_2=setdiff(ind_m_2,ind_t_2) #the patients whose counts in malignant isoform significantly changed

write.csv(ind_2,"part2.csv"); #write into a file
#because of the poor convergence there are much less patients identified here. This is largely because of poor mixing of the MCMC chain
#mixing can be improved by an appropriate choice of the proposal

