/*******************************************************************************
 * Notes:
 * 0) this program is to simulate SNPs with no effect throught generations to
 *    reach drift-mutation equilibrium.
 * 1) mmut is mean number of mutations (mean of a poisson distribution) plus
 *    10 times of the SD.  When lambda is large, poisson distribution approx.
 *    normal distribution.
 * 2) using a comb sub to comb out fixed loci.  This greatly speeded the prog.
 * 3) boundary is checked all the time, but only at the end of splicing or 
 *    insertion, easing coding and gaining speed.
 * 4) SNP output is haplotype majored, i.e., it's haplotype by haplotype
 * 5) Modification to include float mutation rate (2011-.
 *
 *                                       by Xijiang Yu, Dec. 5, 2010
 ******************************************************************************/
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<map>
#include<cstdlib>
#include<algorithm>
#include<ctime>
#include<cmath>
#include<mkl.h>
#include<stdexcept>
#include<sys/stat.h>
using namespace std;

const int max_mutations (1.2e5);
const int tmp_vec_len   (2.4e5);
const int bp_per_morgan (1e8);
const int method_p      (0); //(VSL_RNG_METHOD_POISSON_PTPE); //back compatible
const int method_iu     (0); //(VSL_RNG_METHOD_UNIFORM_STD);
const int comb_period   (1000);
const int timer_period  (500);
const int return_period (timer_period*10);
int       tvec[tmp_vec_len];

int insert(int *a, int m, int *b, int n)
//Insert b into a.  a is of length max_mutations, current length m;
// b of length n holds new mutations
// This sorting program will remove mutations that's same as an existing one.
// Since the uniform range is [0, 1e8], it's not easy to get repeats in new muts
{
  int i, j, k;
  i=j=k=0;

  while(j<n && i<m)
    {
      if(a[i]==b[j])
	{
	  ++j;
	  continue;
	}
      tvec[k++]=(a[i]<b[j])?a[i++]:b[j++];
    }
  while(i<m) tvec[k++]=a[i++];
  while(j<n) tvec[k++]=b[j++];
  if(k>=max_mutations)
    throw runtime_error("Please increase max_mutations in the source code.");
  for(i=0; i<k; ++i) a[i]=tvec[i];
  return k;  //new haplotype length
}

int splice(int seg, int *a, int la, int *b, int lb, int *c, int lc, int *d)
//splice a of length la and  (a, b, c, & d are all sorted in ascending order)
//       b .. ...... lb into
//       d .. ...... ld afterworth, with splicing information
//       c .. ...... lc, ld is returned
//  seg is the initial haplotpye, either 0 for a or 1 for b.
{
  int ia, ib, ic, ld, i;
  ia=ib=ic=ld=0;
  while(ic<lc)
    {
      if(seg)
	{
	  while(ia<la && a[ia]<c[ic]) ++ia;
	  while(ib<lb && b[ib]<c[ic]) tvec[ld++]=b[ib++];
	}
      else
	{
	  while(ia<la && a[ia]<c[ic]) tvec[ld++]=a[ia++];
	  while(ib<lb && b[ib]<c[ic]) ++ib;
	}
      ++ic;
      seg=seg?0:1;
    }
  if(seg) while(ib<lb) tvec[ld++]=b[ib++];
  else    while(ia<la) tvec[ld++]=a[ia++];
  if(ld>=max_mutations)
    throw runtime_error("Please increase max_mutations in the source code.");
  for(i=0; i<ld; ++i) d[i]=tvec[i];
  return ld;  //spliced size
}

//comb out fixed loci
int comb(int *v, int n, int th, map<int, int> &ref)
{
  int i, c[n], m(0);
  // modified 2011.Oct.13, 0 freq loci are also removed
  for(i=0; i<n; ++i) if(ref[v[i]]<th && ref[v[i]]>0) c[m++]=v[i];
  for(i=0; i<m; ++i) v[i]=c[i];
  return m;
}

//output SNPs
void ooput(ostream &oo, int *vec, map<int,int> ref)
{
  map<int,int>::iterator it;
  int i;
  for(i=0, it=ref.begin(); it!=ref.end(); ++it)
    {
      if(vec[i]!=it->first) oo<<0;
      else
	{
	  oo<<1;
	  ++i;
	}
    }
  oo<<'\n';
}

int main(int argc, char*argv[])
{
  if(argc!=5)
    {
      cerr<<"Usage: "<<argv[0]<<" ne ng chr mr"                  <<endl;
      cerr<<endl;
      cerr<<"    ne : integer, effective population size."       <<endl;
      cerr<<"    ng : integer, number of evolving generations."  <<endl;
      cerr<<"    chr: double, chromosome length (in M[organ])."  <<endl;
      cerr<<"    mr : double, mutation rate/meiosis/Morgan."     <<endl;
      cerr<<"    Results will be piped to pgqtl or pdrop"        <<endl;
      cerr<<endl;
      return 1;
    }
  //parameters
  int      ne (atoi(argv[1]));
  int      ng (atoi(argv[2]));
  double   chr(atof(argv[3]));
  double   mr (atof(argv[4]));

  //simulation storage
  int SNPs[2][ne][2][max_mutations], *ptr;
  int parent(0), offspring(1);
  int igrt, iid, ihap, nhap, isnp, status(0);
  int pama[ne], phap[ne]; //phap is the initial parent haplotype
  int nSNP[2][ne][2];         // # mutation on each haplotype
  int poisson[ne];
  int    mxpoint, *point;   // max # and position of mutations/cross-overs
  int    ttbp;              // total # of base pairs
  int    i, k, mv, ml, top, bottom;
  double lambda_m,memory;

  //random number related
  VSLStreamStatePtr stream;
  char   xStream[] = "vsl_Stream_State";
  char   aSeed[]   = "a.seed";
  int    iseed(time(0));
  struct stat  f_info;

  //Different loci index
  map<int,int>           loci;
  map<int,int>::iterator ig;

  //check parameters:
  if(ne<4)
    {
      clog<<"Population: [10, inf].\n";
      ++status;
    }
  if(ng<100)
    {
      clog<<"Generations: [100, integer limit].\n";
      ++status;
    }
  if(chr<0.00001 || chr>4)
    {
      clog<<"Chromosome length: [.00001, 4].\n";
      ++status;
    }
  if(mr<.1 || mr>10)
    {
      clog<<"Mutation rate: [.1, 10].\n";
      ++status;
    }
  if(status) return 2;

  clog<<"\n>>>>> Generating Ideal Population <<<<<\n";
  //Program initialization
  clog<<"\n--> Initializing of parameters & storage ... "<<flush;
  lambda_m = chr*mr;  //parameter for mutation, which follows poisson distr.
  if(mr<1)  mxpoint = int(ne*2*chr      + 10*sqrt(ne*2*chr));
  else      mxpoint = int(ne*2*lambda_m + 10*sqrt(ne*2*lambda_m));
  point=new int[mxpoint];

  ttbp=chr*bp_per_morgan;
  for(i=0; i<2; ++i)
    for(iid=0; iid<ne; ++iid)
      for(ihap=0; ihap<2; ++ihap)
	  nSNP[i][iid][ihap]=0;
  nhap=ne*2;
  clog<<"done."<<endl;

  //Memory usage
  memory  = ne * max_mutations * 2 * 2 * sizeof(int); //mutation storage
  memory += ne * 2 * sizeof(int);                     //sample pa & ma
  memory += ne * 2 * 2 * sizeof(int) * 2 * 2;         //# of mutations
  memory += ne * sizeof(int);                         //# cross-over/muts
  memory += mxpoint * sizeof(int);                    //points storage
  memory += max_mutations * 2 * sizeof(int);     //for insertion and combing
  memory /= (1024*1024);
  memory += 67;                       //Rest memory usage, checked with `top`
  clog<<"--------------> Approximate memory usage ... ";
  clog<<ceil(memory)<<" M\n"<<endl;

  //RNG initialization
  clog<<"--> Initializing random number generator ... "<<flush;
  if(stat(xStream, &f_info)==0)
    {
      vslLoadStreamF(&stream, xStream);
      clog<<"with previous stream.\n"<<endl;
    }
  else{
    if(stat(aSeed, &f_info)==0)
      {
	ifstream fin(aSeed);
	fin>>iseed;
	clog<<"with given integer seed.\n"<<endl;
      }
    vslNewStream(&stream, VSL_BRNG_MT19937, iseed);
    clog<<"with current time as seed.\n"<<endl;
  }

  //Let's create ...
  clog<<"----------------------------> Simulating ..."<<endl;
  //  clog<<"     ";
  for(igrt=0; igrt<ng; ++igrt)
    {
      //create mutations in parents
      for(ihap=0; ihap<2; ++ihap)
	{
	  //how many mutations
	  status = viRngPoisson(method_p, stream, ne, poisson, lambda_m);
	  //where are they
	  status = viRngUniform(method_iu, stream, mxpoint, point, 0, ttbp);
	  //insert the mutations
	  for(iid=i=0; iid<ne; ++iid)
	    {
	      sort(&point[i], &point[i]+poisson[iid]);
	      try
		{
		  nSNP[parent][iid][ihap]
		    =insert(SNPs[parent][iid][ihap],
			    nSNP[parent][iid][ihap],
			    &point[i], poisson[iid]);
		}
	      catch(runtime_error e)
		{
		  cerr<<e.what()<<endl;
		  return 3;
		}
	      i+=poisson[iid];
	    }
	}
      //creating next generation
      //splicing
      for(ihap=0; ihap<2; ++ihap)
	{
	  //sample pa/ma from parent generation
	  top    = (ihap)?ne:(ne/2);
	  bottom = (ihap)?(ne/2):0;
	  status = viRngUniform(method_iu, stream, ne, pama, bottom, top);
	  //the first segment to drop
	  status = viRngUniform(method_iu, stream, ne, phap, 0, 2);
	  //number of cross-overs
	  status = viRngPoisson(method_p , stream, ne, poisson, chr);
	  //where are they
	  status = viRngUniform(method_iu, stream, mxpoint, point, 0, ttbp);
	  for(iid=i=0; iid<ne; ++iid)
	    {
	      sort(&point[i], &point[i]+poisson[iid]);
	      try
		{
		  nSNP[offspring][iid][ihap]
		    =splice(phap[iid],
			    SNPs[parent][pama[iid]][0],
			    nSNP[parent][pama[iid]][0],
			    SNPs[parent][pama[iid]][1],
			    nSNP[parent][pama[iid]][1],
			    &point[i], poisson[iid],
			    SNPs[offspring][iid][ihap]);
		}
	      catch(runtime_error e)
		{
		  cerr<<e.what()<<endl;
		  return 3;
		}
	      i+=poisson[iid];
	    }
	}
      //switch generation
      status    = parent;
      parent    = offspring;
      offspring = status;

      //comb the loci, i.e. remove fixed loci.  500 can be changed to speed up
      if(!((igrt+1)%comb_period)||igrt==ng-1)
	{
	  for(iid=0; iid<ne; ++iid)
	    for(ihap=0; ihap<2; ++ihap)
	      {
		ptr = SNPs[parent][iid][ihap];
		k   = nSNP[parent][iid][ihap];
		for(isnp=0; isnp<k; ++isnp)
		  ++loci[ptr[isnp]];
	      }
	  for(iid=0; iid<ne; ++iid)
	    for(ihap=0; ihap<2; ++ihap)
	      nSNP[parent][iid][ihap] = comb(
					     SNPs[parent][iid][ihap],
					     nSNP[parent][iid][ihap],
					     nhap,loci);
	  loci.clear();
	}
      // program timer
      if((igrt+1)%timer_period  == 0)
	clog<<setw(6)<<igrt+1;
      if((igrt+1)%return_period == 0)
	clog<<endl;
    }
  clog<<endl;

  //post simulation
  for(iid=mv=0; iid<ne; ++iid)
    for(ihap=0; ihap<2; ++ihap)
      {
	ptr = SNPs[parent][iid][ihap];
	k   = nSNP[parent][iid][ihap];
	for(isnp=0; isnp<k; ++isnp)
	  ++loci[ptr[isnp]];
	if(mv<k) mv=k;
      }
  clog<<"---> Max haplotype storage vector length ... "<<mv<<endl;

  //before piping to gdrop or other program, write the stream file.
  vslSaveStreamF(stream, xStream);

  //parameters to be piped to next program
  cout<<ne<<endl<<loci.size()<<endl;
  // modified <2011.Oct.13>, output the genotypes first
  //create genotype file
  for(iid=0; iid<ne; ++iid)
    for(ihap=0; ihap<2; ++ihap)
      ooput(cout, SNPs[parent][iid][ihap], loci);

  //creating linkage map
  for(ml=0, ig=loci.begin(); ig!=loci.end(); ++ig)
    if(ig->second<nhap)
      cout<<setw(9)<<ig->first<<setw(5)<<ig->second<<'\n';

  if(point) delete []point;
  clog<<"----> Finalizing random number generator ..."<<flush;
  vslDeleteStream(&stream);
  clog<<" done.\n"<<endl;

  return 0;
}
