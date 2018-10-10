#include <mpi.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>


using namespace std;
double f(double x) {
  double y;
  int i;
  y = x;
  for(i = 1;i<=10;i++){
    y += sin(x*i)/pow(2,i);
  }
  return y;
}

int main(int argc, char **argv){
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);					//the edge size of the array A
    int p;									//number of processors
    int id;									//processor id

    MPI_Init(&argc, &argv);					//initialize MPI


    MPI_Comm_size(MPI_COMM_WORLD,&p);		//the total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD,&id);		//processor id

    int process_rown;
    int process_coln;
    if (p==1) {
        process_rown = 1;
        process_coln = 1;
    }
    else{
        process_rown = sqrt(round(double(m)/double(n)*double(p)));
        process_coln = sqrt(round(double(n)/double(m)*double(p)));
    }



    int ms = (m%process_rown==0)?(m/process_rown):floor(m/process_rown)+1;			//the size of every subarray on each processor
    int ns = (n%process_coln==0)?(n/process_coln):floor(n/process_coln)+1;

    double **F_map = new double*[ms+2];
    for(int i = 0; i < ms+2;i++){
        F_map[i] = new double[ns+2];
    }

    for(int i = 0;i<ms+2;i++){
        for(int j = 0;j<ns+2;j++){
            F_map[i][j] = 0;
        }
    }


    double **A_local = new double*[ms+2];
    for(int i = 0; i < ms+2;i++){
        A_local[i] = new double[ns+2];
    }

    for(int i = 0;i<ms+2;i++){
        for(int j = 0;j<ns+2;j++){
            A_local[i][j] = 0;
        }
    }

    double **z = new double*[ms];
    for(int i = 0; i < ms;i++){
        z[i] = new double[ns];
    }

    for(int i = 0;i<ms;i++){
        for(int j = 0;j<ns;j++){
            z[i][j] = 0;
        }
    }
    

    int pid_x;								//the row id of each processor
    int pid_y;								//the column id of each processor

    double starttime;						//the start time of 10 times iterations
    double endtime;							//the end time of 10 times iterations

    double endtime1;
    double endtime2;

    double elapsedtime;						//the elapsed time of this MPI program

    double elapsedtime1;
    double elapsedtime2;

    int counts=0;							//count the number of the iterations

    MPI_Status status;						//MPI status for MPI_Recv function

    pid_x=floor(id/process_coln);
    pid_y=id-process_coln*pid_x;
    //cout<<"id"<<id<<"pid_x"<<pid_x<<"pid_y"<<pid_y<<endl;
    int ireal;
    int jreal;
    bool UP = 1;
    bool DW = 1;
    bool LF = 1;
    bool RT = 1;

    int i_length;
    int j_length;
    int i_start;
    int j_start;

    double send_to_top[ns];             
    double send_to_left[ms];             
    double send_to_right[ms]; 
    double send_to_bottom[ns]; 

    double receive_from_bottom[ns];       //receive a one-line array from bottom subarray
    double receive_from_top[ns];        //receive a one-line array to right subarray
    double receive_from_right[ms];      //rceive a single value to bottomright subarray
    double receive_from_left[ms];

    //initialize array A
    if(pid_x!=process_rown-1&&pid_y!=process_coln-1){
        for(int i=1;i<ms+1;i++){
            for(int j=1;j<ns+1;j++){
                ireal = pid_x*ms+i-1;
                jreal = pid_y*ns+j-1;
                A_local[i][j]=ireal*sin(ireal)+jreal*cos(jreal)+sqrt(ireal+jreal);
            }
        }
    }
    else if(pid_x==process_rown-1&&pid_y!=process_coln-1){
        for(int i=1;i<m+ms-ms*process_rown+1;i++){
            for(int j=1;j<ms+1;j++){
                ireal = pid_x*ms+i-1;
                jreal = pid_y*ns+j-1;
                A_local[i][j]=ireal*sin(ireal)+jreal*cos(jreal)+sqrt(ireal+jreal);
            }
        }
    }
    else if(pid_x!=process_rown-1&&pid_y==process_coln-1){
        for(int i=1;i<ms+1;i++){
            for(int j=1;j<n+ns-ns*process_coln+1;j++){
                ireal = pid_x*ms+i-1;
                jreal = pid_y*ns+j-1;
                A_local[i][j]=ireal*sin(ireal)+jreal*cos(jreal)+sqrt(ireal+jreal);
            }
        }
    }
    else{
        for(int i=1;i<m+ms-ms*process_rown+1;i++){
            for(int j=1;j<n+ns-ns*process_coln+1;j++){
                ireal = pid_x*ms+i-1;
                jreal = pid_y*ns+j-1;
                A_local[i][j]=ireal*sin(ireal)+jreal*cos(jreal)+sqrt(ireal+jreal);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int id_snip = 0;
    int i_length_in = i_length;
    int j_length_in = j_length;
    if(id==0){
        starttime=MPI_Wtime();
        cout<<"The 10 tims iterations start at "<<starttime<<endl;
    }



    
    //10 times iterations
    while(counts<10){
        i_length = ms;
        j_length = ns;
        i_start = 1;
        j_start = 1;
        //cout<<"send imformation to the top and left 3 subarrays"<<id<<"("<<pid_x<<","<<pid_y<<")"<<endl;
        //send imformation to the top and left 3 subarrays
        counts++;
        //if(id == id_snip) cout<<"id"<<id<<"round"<<counts<<endl;
        if (pid_x==0){
            UP = 0;
            i_start = 2;
        }
        if (pid_y==0){
            LF = 0;
            j_start = 2;
        }   

        if (pid_x==process_rown-1){
            DW = 0;
            i_length = m+ms-ms*process_rown;
        }
        if (pid_y==process_coln-1){
            RT = 0;
            j_length = n+ns-ns*process_coln;
        }   


        //if(id == id_snip) cout<<"id"<<id<<"LF"<<LF<<"RT"<<RT<<"UP"<<UP<<"DW"<<DW<<endl;
//////////////////SEND Left Rec Right//////////
        if (LF){
            for(int i=1;i<i_length+1;i++){
                send_to_left[i-1]=A_local[i][1];
            }
            MPI_Ssend(&send_to_left,i_length,MPI_DOUBLE,pid_x*process_coln+(pid_y-1),0,MPI_COMM_WORLD);
        }
        
        if(RT){
            MPI_Recv(&receive_from_right,i_length,MPI_DOUBLE,pid_x*process_coln+(pid_y+1),0,MPI_COMM_WORLD,&status);
//////////////////SEND Right Rec Left//////////
            for(int i=1;i<i_length+1;i++){
                send_to_right[i-1]=A_local[i][ns];
            }
            MPI_Ssend(&send_to_right,i_length,MPI_DOUBLE,pid_x*process_coln+(pid_y+1),1,MPI_COMM_WORLD);
        }

        if(LF){
            MPI_Recv(&receive_from_left,i_length,MPI_DOUBLE,pid_x*process_coln+(pid_y-1),1,MPI_COMM_WORLD,&status);
        }
//////////////////SEND UP Rec DOWN//////////
        if (UP){
            for(int j=1;j<j_length+1;j++){
                send_to_top[j-1]=A_local[1][j];
            }
            MPI_Ssend(&send_to_top,j_length,MPI_DOUBLE,(pid_x-1)*process_coln+pid_y,2,MPI_COMM_WORLD);

        }

        if(DW){
            MPI_Recv(&receive_from_bottom,j_length,MPI_DOUBLE,(pid_x+1)*process_coln+pid_y,2,MPI_COMM_WORLD,&status);
//////////////////SEND DOWN Rec UP//////////
            for(int j=1;j<j_length+1;j++){
                send_to_bottom[j-1]=A_local[ms][j];
            }
            MPI_Ssend(&send_to_bottom,j_length,MPI_DOUBLE,(pid_x+1)*process_coln+pid_y,3,MPI_COMM_WORLD);
        }
        
        if(UP){
            MPI_Recv(&receive_from_top,j_length,MPI_DOUBLE,(pid_x-1)*process_coln+pid_y,3,MPI_COMM_WORLD,&status);
        }
        i_length_in = i_length;
        j_length_in = j_length;
        if(RT)  j_length_in = j_length+1; //Now ns+1 or n+ns-ns*process_coln
        if(DW)  i_length_in = i_length+1; //Now ms+1 or m+ms-ms*process_rown

        for(int i = i_start;i<i_length_in;i++){
            if(LF) A_local[i][0] = receive_from_left[i-1];
            if(RT) A_local[i][ns+1] = receive_from_right[i-1];
        }
        for(int j = j_start;j<j_length_in;j++){
            if(UP) A_local[0][j] = receive_from_top[j-1];
            if(DW) A_local[ms+1][j] = receive_from_bottom[j-1];
        }




        //endtime1 = MPI_Wtime();
        for(int i = i_start-1;i<i_length_in+1;i++){
            for(int j = j_start-1;j<j_length_in+1;j++){
                F_map[i][j] = f(A_local[i][j]);
            }
        }
        //endtime2=MPI_Wtime();
        for(int i = i_start;i<i_length_in;i++){
            for(int j = j_start;j<j_length_in;j++){
                double temp = (F_map[i+1][j]+F_map[i][j+1]+F_map[i][j]+F_map[i-1][j]+F_map[i][j-1])/5;
                z[i-1][j-1] = temp>100?100:(temp<-100?-100:temp);
            }
        }
        
        for(int i = i_start;i<i_length_in;i++){
            for(int j = j_start;j<j_length_in;j++){
                A_local[i][j] = z[i-1][j-1];
            }
        }

    }
    
    double sum_local;
    double sum_local_square;
    double sum_temp;
    double sum_temp_square;
    double sum;
    double sum_square;

    // get the sum of A;
    // fisrt get the local sum of each processor
    for(int i=1;i<i_length+1;i++){
        for(int j=1;j<j_length+1;j++){
            sum_local=sum_local+A_local[i][j];
            sum_local_square = sum_local_square+A_local[i][j]*A_local[i][j];
        }
    }

    //send local sum to root processor and get the sum of array A
    if(pid_y!=0){
        MPI_Ssend(&sum_local,1,MPI_DOUBLE,pid_x*process_coln,pid_y,MPI_COMM_WORLD);
        MPI_Ssend(&sum_local_square,1,MPI_DOUBLE,pid_x*process_coln,pid_y+p,MPI_COMM_WORLD);
    }
    if(pid_y==0){
        for(int j=1;j<process_coln;j++){
            MPI_Recv(&sum_temp,1,MPI_DOUBLE,pid_x*process_coln+j,j,MPI_COMM_WORLD,&status);
            MPI_Recv(&sum_temp_square,1,MPI_DOUBLE,pid_x*process_coln+j,j+p,MPI_COMM_WORLD,&status);
            sum_local=sum_local+sum_temp;
            sum_local_square = sum_local_square+sum_temp_square;
        }
    }

    if(pid_y==0&&pid_x!=0){
        MPI_Ssend(&sum_local,1,MPI_DOUBLE,0,pid_x,MPI_COMM_WORLD);
        MPI_Ssend(&sum_local_square,1,MPI_DOUBLE,0,pid_x+2*p,MPI_COMM_WORLD);
    }
    if(pid_y==0&&pid_x==0){
        for(int i=1;i<process_rown;i++){
            MPI_Recv(&sum_temp,1,MPI_DOUBLE,i*process_coln,i,MPI_COMM_WORLD,&status);
            MPI_Recv(&sum_temp_square,1,MPI_DOUBLE,i*process_coln,i+2*p,MPI_COMM_WORLD,&status);
            sum_local=sum_local+sum_temp;
            sum_local_square = sum_local_square+sum_temp_square;
        }
        sum = sum_local;
        sum_square = sum_local_square;
        endtime=MPI_Wtime();
        elapsedtime=endtime-starttime;
        //elapsedtime1=endtime1-starttime;
        //elapsedtime2=endtime2-endtime1;
        cout<<"The elapsed time is  "<<elapsedtime<<endl;
        //cout<<"The endtime1-starttime is  "<<elapsedtime1<<endl;
        //cout<<"The endtime2-endtime1 is  "<<elapsedtime2<<endl;


        cout<<"The sum of array A is  "<<sum<<endl;
        cout<<"The sum square of array A is  "<<sum_square<<endl;
    }
/*
    if(pid_y!=0){
        MPI_Ssend(&sum_local_square,1,MPI_DOUBLE,pid_x*process_coln,pid_y+p,MPI_COMM_WORLD);
    }
    if(pid_y==0){
        for(int j=1;j<process_coln;j++){
            MPI_Recv(&sum_temp_square,1,MPI_DOUBLE,pid_x*process_coln+j,j+p,MPI_COMM_WORLD,&status);
            sum_local=sum_local+sum_temp;
            sum_local_square = sum_local_square+sum_temp_square;
        }
    }

    if(pid_y==0&&pid_x!=0){
        MPI_Ssend(&sum_local_square,1,MPI_DOUBLE,0,pid_x+2*p,MPI_COMM_WORLD);
    }
    if(pid_y==0&&pid_x==0){
        for(int i=1;i<process_rown;i++){
            MPI_Recv(&sum_temp_square,1,MPI_DOUBLE,i*process_coln,i+2*p,MPI_COMM_WORLD,&status);
            sum_local_square = sum_local_square+sum_temp_square;
        }
        sum_square = sum_local_square;
        endtime=MPI_Wtime();
        elapsedtime=endtime-starttime;
        cout<<"The elapsed time is  "<<elapsedtime<<endl;
        cout<<"The sum of array A is  "<<sum<<endl;
    }
*/
    //cout<<"send verification value center of array to root processor"<<endl;
    //send verification value center of array to root processor

    MPI_Finalize();

    return 0;
}

