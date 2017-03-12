#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "omp.h"
struct bigNum
{
	int *num;
	int len;
};

struct bigNum input(FILE* fp,int size)
{
	struct bigNum b;	
	b.num=malloc(sizeof(int)*(size));
	char ch;
	int l=0;
	while((ch=fgetc(fp))!=EOF && ch!='\n')
	{
		if(l==(size))
		{
			b.num=realloc(b.num,sizeof(int)*((size)+=4));
			if(!b.num)//if integer pointer couldn't be reallocated
				return b;
		}
		b.num[l++]=ch-'0';
	}
	b.len=l;
	return b;
}

void print(struct bigNum a)
{
	int i=a.len-1;
	while(a.num[i]==0&&i>0)
		i--;
	for(;i>=0;i--)
		printf("%d",a.num[i]);
	printf("\n");
}

int exp2k(int m,int n)
{
	int p=1;
	while(p<=m || p<=n)
		p*=2;
	return p;
}

void digitize(struct bigNum *a)
{
	int i,n=(*a).len,carry=0;
	for(i=0;i<n;i++)
	{
		(*a).num[i]+=carry;		
		if((*a).num[i]<0)
			carry=(-((*a).num[i]+1)/10+1);
		else
			carry=(*a).num[i]/10;
		(*a).num[i]=(*a).num[i]-carry*10;
	}
}

void karatsuba(struct bigNum *a,struct bigNum *b,struct bigNum *c,int pow_2k)
{
	if(pow_2k==1)
	{
		(*c).num[0]=(*c).num[1]=0;
		(*c).num[0]=(*a).num[0]*(*b).num[0];
	}
	else
	{
		int m2=pow_2k/2,i;
		struct bigNum a1,a2,b1,b2,p1,p2,p3,sa,sb;
		a1.len=a2.len=b1.len=b2.len=sa.len=sb.len=p1.len=p2.len=p3.len=m2;
		a1.num=&((*a).num[0]);//Lower bits a
		a2.num=&((*a).num[m2]);//Higher bits a
		b1.num=&((*b).num[0]);//Lower bits b
		b2.num=&((*b).num[m2]);//Higher bits b
		//Using extra space in (*c).num to store the sum of 2 halves
		sa.num=&((*c).num[5*pow_2k]);//Sum of halves of a
		sb.num=&((*c).num[5*pow_2k+m2]);//Sum of halves of b		
		for(i=0;i<m2;i++)
		{
			sa.num[i]=a1.num[i]+a2.num[i];
			sb.num[i]=b1.num[i]+b2.num[i];
		}
		p1.num=&((*c).num[pow_2k]);//p1*10^(2^k)
		p3.num=&((*c).num[2*pow_2k]);//p3 stored using extra space to avoid segmentation faults for large arrays
		p2.num=&((*c).num[0]);//p2*10^0
		karatsuba(&a1,&b1,&p2,m2);
		karatsuba(&a2,&b2,&p1,m2);
		karatsuba(&sa,&sb,&p3,m2);
		//#pragma omp parallel for
		for(i=0;i<pow_2k;i++)
			p3.num[i]=p3.num[i]-p2.num[i]-p1.num[i];//p3-p2-p1
		//#pragma omp parallel for
		for(i=0;i<pow_2k;i++)
			(*c).num[i+m2]+=p3.num[i];//(p3-p2-p1)*10^((2^k)/2)
	}
}

double exec_time(struct timespec start,struct timespec end)
{
	double d=end.tv_sec*1000000000+end.tv_nsec;
	d=d-(start.tv_sec*1000000000+start.tv_nsec);
	return d/1000000;
}

int main()
{
	struct bigNum a=input(stdin,4);
	struct bigNum b=input(stdin,4);
	struct timespec start,end;
	int pow_2k=exp2k(a.len,b.len);//smallest power of 2 greater than both a.len & b.len
	//preformatting inputs
	int len=pow_2k,i;
	struct bigNum x,y;//x&y stores reverse of a&b
	x.num=malloc(sizeof(int)*len);
	y.num=malloc(sizeof(int)*len);
	x.len=y.len=len;
	#pragma omp parallel for // uncomment this to view execution without parallelism
	for(i=0;i<len;i++)//Reverse the array and store, 0 index represesnts units place of number
	{
		if(i<a.len)
			x.num[i]=a.num[a.len-1-i];
		else
			x.num[i]=0;
		if(i<b.len)
			y.num[i]=b.num[b.len-1-i];
		else
			y.num[i]=0;	
	}
	free(a.num);
	free(b.num);
	struct bigNum z;
	z.num=malloc(sizeof(int)*(6*len));//0-2:Actual Product,2-4 extra space to compute p3-p2-p1,5-6 sum of two halves of 2 numbers
	z.len=2*len;
	clock_gettime(CLOCK_REALTIME,&start);
	karatsuba(&x,&y,&z,pow_2k);
	digitize(&z);
	clock_gettime(CLOCK_REALTIME,&end);
	//print(z);not printing output on terminal
	printf("Execution Time:%lfms\n",exec_time(start,end));
	free(x.num);
	free(y.num);
	free(z.num);
	return 0;
}
