#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include <deque>
#include <queue>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include "f.h" 

using namespace std;

struct interval{
	double a;
	double b;
	double fa;
	double fb;
	double maximum;
};

struct comparator{
    bool operator()(const interval &a,const interval &b){
        return a.maximum<b.maximum;
    }
};


int main(int argc, char **argv){
	
	double a=stod(argv[1]);	
	double b=stod(argv[2]);	
	double e=stod(argv[3]);	
	double s=stod(argv[4]);	
	int n;									
	int id;									
	double start_time;						 
	double end_time;						
	double total_time;					    

	priority_queue<interval,vector<interval>,compare> MAXs;
	
	double Maxin;
    vector<bool>done;
    bool finish=false;

	start_time = omp_get_wtime();

	#pragma omp parallel private(id)	
	{
		n = omp_get_num_threads();
    	id = omp_get_thread_num();		
    	#pragma omp single
    	{
    		interval sing_interval;
    		sing_interval.a=a;
			sing_interval.b=b;
			sing_interval.fa=f(a);
			sing_interval.fb=f(b);
			sing_interval.maximum=(sing_interval.fa+sing_interval.fb+s*(b-a))/2;
			done.resize(n);
			for(int i=0;i<n;i++){
				done[i]=false;
			}

			MAXs.push(sing_interval);
			Maxin=max(sing_interval.fa,sing_interval.fb);
    	}
    	
    	if(id==0){				
			//Master
			interval local_interval;
			vector<bool> master_done;
			double local_Maxin;
    		interval local_interval1;
    		interval local_interval2;
    		bool local_first=false;

    		while(!local_first){
    			#pragma omp critical
    			{
    				if(!MAXs.empty()){
    					local_interval=MAXs.top();	
    					MAXs.pop();
    					local_first=true;
    				}
    			}    			
    		}

			while(1){
				bool local_flag=true;
				bool local_update=false;				

    			if(Maxin+e>local_interval.maximum){
    				for(size_t i=1;i<done.size();i++){
    					local_flag=local_flag&master_done[i];
    				}
    				if(local_flag){
    					//omp_set_lock(&lock);
    					#pragma omp critical
    					{
    						finish=true;
    					}    					
    					//omp_unset_lock(&lock);
    					break;
    				}
    			}
    			else{
    				local_Maxin=max(max(local_interval.fa,local_interval.fb),Maxin);

					local_update=true;
					local_interval1.a=local_interval.a;
					local_interval1.b=(local_interval.a+local_interval.b)/2;
					local_interval1.fa=local_interval.fa;
					local_interval1.fb=f(local_interval1.b);
					local_interval1.maximum=(local_interval1.fa+local_interval1.fb+s*(local_interval1.b-local_interval1.a))/2;//(fa+fb+s*(b-a))/2
					local_interval2.a=(local_interval.a+local_interval.b)/2;
					local_interval2.b=local_interval.b;
					local_interval2.fa=local_interval1.fb;
					local_interval2.fb=local_interval.fb;
					local_interval2.maximum=(local_interval2.fa+local_interval2.fb+s*(local_interval2.b-local_interval2.a))/2				
				}
				if(local_update){
					#pragma omp critical
    				{
    					MAXs.push(local_interval1);
						MAXs.push(local_interval2);
						local_interval=MAXs.top();	
    					MAXs.pop();
    					master_done=done;
    					Maxin=local_Maxin;					
    				}
    				
    			}
    			else{
    				#pragma omp critical
    				{
						local_interval=MAXs.top();	
    					MAXs.pop();
    					master_done=done;					
    				}
    			}    								

    			//count++;		
				//printf("Finishing one time iteration %d\n",count);
    		}		
    	}
    	else{
    		
    		interval local_interval;
    		double local_Maxin;
    		interval local_interval1;
    		interval local_interval2;
    		bool local_done;
    		bool local_first=false;

    		while(!local_first){
    			#pragma omp critical
    			{    			
    				if(!MAXs.empty()){
    					local_interval=MAXs.top();	
    					MAXs.pop();
    					local_first=true;
    				}
    			}    			
    		}
    		//printf("The total thread number is  %d.The thread %d starts\n",n ,id);
    		
    		while(!finish){  
    			bool local_update=false;			

   				if(Maxin+e<local_interval.maximum){			
    				local_done=false;			
					
					local_Maxin=max(max(local_interval.fa,local_interval.fb),Maxin);

					local_update=true;
					local_interval1.a=local_interval.a;
					local_interval1.b=(local_interval.a+local_interval.b)/2;
					local_interval1.fa=local_interval.fa;
					local_interval1.fb=f(local_interval1.b);
					local_interval1.maximum=(local_interval1.fa+local_interval1.fb+s*(local_interval1.b-local_interval1.a))/2;
					
                    local_interval2.a=(local_interval.a+local_interval.b)/2;
					local_interval2.b=local_interval.b;
					local_interval2.fa=local_interval1.fb;
					local_interval2.fb=local_interval.fb;
					local_interval2.maximum=(local_interval2.fa+local_interval2.fb+s*(local_interval2.b-local_interval2.a))/2;    					

    			}
    			else{    					
    				local_done=true;
    			}
    			if(!local_update){
    				#pragma omp critical
    				{    					
    					local_interval=MAXs.top();	
    					MAXs.pop();
    					done[id]=local_done;
    				}    									
    			}
    			else{
    				#pragma omp critical
    				{
    					MAXs.push(local_interval1);
						MAXs.push(local_interval2);
						Maxin=local_Maxin;
						done[id]=local_done;
						local_interval=MAXs.top();	
    					MAXs.pop();						
    				}				
				}    							
    		}
    	}

    }
    end_time=omp_get_wtime();
    total_time =end_time-start_time;

    printf("Time:  %f\n",total_time);
	printf("Maximum: %f\n",M);
	return 0;
}
