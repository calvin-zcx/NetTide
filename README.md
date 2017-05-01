# NetTide
@inproceedings{Zang2016BeyondST,
  title={Beyond Sigmoids: The NetTide Model for Social Network Growth, and Its Applications},
  author={Chengxi Zang and Peng Cui and Christos Faloutsos},
  booktitle={KDD},
  year={2016}
}

More details in the upcoming journal version.

More info on the Social Dynamics and beyond: http://media.cs.tsinghua.edu.cn/~multimedia/cuipeng/

Author: Chengxi Zang
Date: 2017-05-01

C++: 
Goal: Generate realistic growth dynamics, both for node and link, of social networks.
The growth dynamics of node and link are captured by NETTIDE equations:\n
	*Node adoption:
		n-dot = beta / t^theta * n * (N - n )
	*Link build between infected-infected and infected-newly infected:
		e-dot = beta_prime / t^theta * n * ( alpha * (n-1)^gamma - e/n ) + 2*n-dot
	*Assumption:
		beta and beta_prime : growth rate of nodes and links;
		1/t^theta : fizzling effect
		the alpha and gamma describe the average number of friends in the ego network.
Two generators:
	*NetTide-Process: captures the micro-level stochastic interactions within a network.
	*NetTide-Survival: captures the growth dynamics from hazard rates (hazard process).
Input 
	*Input network, reads the edge file 
	    Data format: blank-separated pairs , e.g 0 1 represents 0->1
		Network: e.g.,
			1.Complete network, random network
			2.Kronecker network
			3.Or any other network
        Nodes have been numbered 0-$n-1$
	*Parameters: beta, theta, beta_prime, alpha, gamma, N
Output:
  *n(t) over time 
	*e(t) over time. 
	*For NetTide-Process, detailed evovling process of nodes and links are recorded.

Matlab:
Goal: Fit the growth dynamics by NetTide model, plot, and regenerate them.
