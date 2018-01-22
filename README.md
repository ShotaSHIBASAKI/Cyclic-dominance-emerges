# Shibasaki & Shimada 2018

The codes are used for the simulations in :
* Shibasaki and Shimada, (2018), Cyclic dominance emerges from 
the two cooperative behaviros in the social amoeba.

# Description  
## Asexual Continuous Replicator.py  
This code is for the model excluding the mating types and the spacial strucutre. 
Here, the continuous replicator dynamics is applied.   
We used odeint function here for simplicity.  
This code shows you that the evolutionary dynamics converge to the interior point  
if the cyclic dominance appears.  
You can find a figure simmilar to Fig. 2.  

## Asexual Discrete Condition.py  
This code check the conditions where the interior equilobrium G is linearly stable
in the asexual model with discrete replicator dynamics.  
If all absolute values of the eigenvalues of linearized matrix L is smaller than 1,  
the interor point G is stable.  
We can compare the results with those in the continuous replicator dynamics.  
You can find a figure simmilar to Fig. S1 a.  

## Asexual Discrete Replicator.py
This code is for the model excluding the mating types and the spacial strucutre. 
Here, the discrete replicator dynamics is applied. 
This code also check whether the the interipr point G is stable or not in two ways.     
This code shows you that the evolutionary dynamics converge to the interior point  
if the cyclic dominance appears.  
You can find a figure simmilar to Fig. S1 b. 

## Sexual Continuous Replicator.py  
This code is for the model including the mating types but without the spacial strucutre. 
Here, the dcontinuous replicator dynamics is applied.     
You can find figures simmilar to Fig. S2.  

## Experimental Data Analysis.py  
This code is for analyzing the experimental results.  
Then, classify the results by using the linear support vector machine  
named LinearSVC in scikit-learn.  
The parameter for LinearSVC is fitted by the grid research.  
You can find figures simillar to Fig. S7.

## Dynamic Environment Sexual Discrete Replicator.py  
This code simulates the evolutionary dynamics of the sexual model  
under the dynamic environement.  
The results of the evolutionary dynamics will be changed according to the environemnt.  
You can find figures similar to Fig. S8.

## Agent-base Sexual Mutation-.c  
This code is for simulating tye agent-based model  
including the mating types and the spacial structure, but excluding the mutaiton.  
Moore neighborhood and the stochastic update process are applied.
You can generate the data as shown in Fig. 4.  
In addition, you can also generate the data as in Fig. S3-6 by arranging this code.

## Agent-base Sexual Mutation+.c  
This code is for simulating tye agent-based model  
including the mating types, the spacial structure, and the mutaiton.  
Mutation occurs in each 500 time steps.  
Moore neighborhood and the stochastic update process are applied.
You can generate the data as shown in Fig. 5.  
