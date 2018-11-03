#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include <deque>
#include <queue>
#include <vector>
#include <stdio.h>
#include <omp.h>

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


double f(double x);

void split_interval(interval &local_interval1,interval &local_interval2, interval local_interval,double s);


int main(int argc, char **argv){
	
	double a=atof(argv[1]);	
	double b=atof(argv[2]);	
	double e=atof(argv[3]);	
	double s=atof(argv[4]);	
	int n;									
	int id;														    
	priority_queue<interval,vector<interval>,comparator> MAXs;
	
	double Maxout;
    vector<bool>done;
    bool finish=false;


    double start_time;                       
    double end_time;                        
    double total_time;
	start_time = omp_get_wtime();

	#pragma omp parallel private(id)	
	{
        //initialization
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
			Maxout=max(sing_interval.fa,sing_interval.fb);
    	}

    	//Master
    	if(id==0){				
			interval local_interval;
			vector<bool> master_done;
			double local_Maxout;
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

    			if(Maxout+e>local_interval.maximum){
    				for(size_t i=1;i<done.size();i++){
    					local_flag=local_flag&master_done[i];
    				}
    				if(local_flag){
    					#pragma omp critical
    					{
    						finish=true;
    					}    					
    					break;
    				}
    			}
    			else{
    				local_Maxout=max(max(local_interval.fa,local_interval.fb),Maxout);
					local_update=true;
					split_interval(local_interval1,local_interval2,local_interval,s);			
				}
				if(local_update){
					#pragma omp critical
    				{
    					MAXs.push(local_interval1);
						MAXs.push(local_interval2);
						local_interval=MAXs.top();	
    					MAXs.pop();
    					master_done=done;
    					Maxout=local_Maxout;					
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
    		}		
    	}
    	else{
    		interval local_interval;
    		double local_Maxout;
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
    		
    		while(!finish){  
    			bool local_update=false;			

   				if(Maxout+e<local_interval.maximum){			
    				local_done=false;			
					local_Maxout=max(max(local_interval.fa,local_interval.fb),Maxout);
					local_update=true;
					split_interval(local_interval1,local_interval2,local_interval,s);
    			}
    			else local_done=true;
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
						Maxout=local_Maxout;
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
	printf("Maximum: %f\n",Maxout);
	return 0;
}


double f(double x) {
    double ans = 0;
    for (int i = 100; i >= 1; --i) {
        double temp = 0;
        for (int j = i; j >= 1; --j) temp += pow(x + 0.5*j, -3.3);
        ans += sin(x + temp) / pow(1.3, i);
    }
    return ans;
}


void split_interval(interval &local_interval1,interval &local_interval2, interval local_interval,double s) {
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
