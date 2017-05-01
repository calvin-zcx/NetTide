# NetTide
Ref:   
@inproceedings{Zang2016BeyondST,<br>
  title={Beyond Sigmoids: The NetTide Model for Social Network Growth, and Its Applications},<br>
 author={Chengxi Zang and Peng Cui and Christos Faloutsos},<br>
  booktitle={KDD},<br>
 year={2016}<br>
}<br>

More details in the upcoming journal version.<br>
More info on the Social Dynamics and beyond: http://media.cs.tsinghua.edu.cn/~multimedia/cuipeng/  <br>

Author: Chengxi Zang<br>
Date: 2017-05-01<br>

C++:   

>Goal: Generate realistic growth dynamics, both for node and link, of social networks.<br>   
>The growth dynamics of node and link are captured by NETTIDE equations:<br>
>>	Node adoption:<br>
>>>		n-dot = beta / t^theta * n * (N - n) 
>>	Link build between infected-infected and infected-newly infected:
>>>		e-dot = beta_prime / t^theta * n * ( alpha * (n-1)^gamma - e/n ) + 2*n-dot 
>>	Assumption: <br>
>>>		beta and beta_prime : growth rate of nodes and links; 
>>>		1/t^theta : fizzling effect 
>>>		the alpha and gamma describe the average number of friends in the ego network.
>Two generators:<br>
>>	NetTide-Process: captures the micro-level stochastic interactions within a network.     
>>	NetTide-Survival: captures the growth dynamics from hazard rates (hazard process).    

>Input:<br>   
>>	Input network, reads the edge file 
>>  Data format: blank-separated pairs , e.g 0 1 represents 0->1 
>>	Network: e.g., 
>>>			1.Complete network, ra ndom network 
>>>			2.Kronecker network 
>>>			3.Or any other network 
>>  Nodes have been numbered 0-$n-1$ 
>>	Parameters: beta, theta, beta_prime, alpha, gamma, N   

>Output: <br>  
>>  n(t) over time  
>>	e(t) over time. 
>>	For NetTide-Process, detailed evovling process of nodes and links are recorded. 

Matlab:      
>Goal: Fit the growth dynamics by NetTide model, plot, and regenerate them. <br>

