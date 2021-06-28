#include "parseVCF.h"
#include "sort.h"



int verbose=0;

#define COMMAND_CONV 0
#define COMMAND_ANALYSIS 1
#define COMMAND_STATS 2
#define COMMAND_FILTER 3

#define COMMAND_CONV_GEN 0
#define COMMAND_CONV_DIP 1
#define COMMAND_CONV_DOSE 2
#define COMMAND_CONV_GL 3
#define COMMAND_CONV_AS 4
#define COMMAND_CONV_SUN 5
#define COMMAND_CONV_AP 6
#define COMMAND_CONV_HAPMAP 7
#define COMMAND_CONV_EIGEN 8
#define COMMAND_CONV_IMPUTE 9
#define COMMAND_CONV_DHETR 10

#define COMMAND_ANALYSIS_AS 0
#define COMMAND_ANALYSIS_GENCALL 1
#define COMMAND_ANALYSIS_BF 2

#define COMMAND_STATS_GENERAL 0
#define COMMAND_STATS_RSQ 1
#define COMMAND_STATS_GENCOUNT 2

void printDouble(double* v, long n){
	long i;
        for(i=0; i<n; i++){fprintf(stderr, "%lf ", v[i]);}
        fprintf(stderr, "\n");
}
void printLong(long* v, long n){
	long i;
	for(i=0; i<n; i++){fprintf(stderr, "%ld ", v[i]);}
	fprintf(stderr, "\n");
}
int isExon(long pos, long* starts, long* ends, long nexon){
	int i;
	if(nexon==0){return 1;}
	for(i=0; i<nexon; i++){
		if(starts[i]<=pos && pos<=ends[i]){
			return 1;
		}
	}
	return 0;
}

long countFields(char* posStr){
    long i;
    long n=0;
    if(strcmp(posStr,"NA")==0){return 0;}
    for(i=0; i<strlen(posStr); i++){
        if(posStr[i]==','){
            n++;
        }
    }
    return n;
}
long splitCSV(char* posStr, long n, long* pos){
    long i;
    long offset=0;
    char posChar[100];
    for(i=0; i<n; i++){
        sscanf(posStr+offset, "%ld,", pos+i);
        sprintf(posChar, "%ld", pos[i]);
        offset += strlen(posChar)+1;
    }
    return 0;
}


void clearLong(long* v, long N){
	long i;
	for(i=0; i<N; i++){v[i] = 0;}
}
long sum(double* x, int n){
	int i;
	long tot=0;
	for(i=0; i<n; i++){
		tot+=x[i];
	}
	return tot;
}

long sumLong(long* x, int n){
        int i;
        long tot=0;
        for(i=0; i<n; i++){
                tot+=x[i];
        }
        return tot;
}

int covCheck(long* x, int n, int th){
	int i;
	for(i=0; i<n; i++){
		if(x[i*2]+x[i*2+1]>th){return 1;}	
	}
	return 0;
}

long sumL(long* x, int n){
	int i;
	long res=0;
	for(i=0; i<n; i++){
		res += x[i];
	}
	return res;
}

void print(double* x, int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<2;j++){
			printf("%lf\t",x[i*3+j]);
		}
		printf("\n");
	}
}

int main(int argc, char** argv){
	
	int m=0, M=0; // another file for Rsq calc
	int* Dip; // dips for Rsq calc
	char** rssAdd;
	
	if(init()==0){fprintf(stderr, "error on parseVCF\n"); return -1;};
        int i;
        if(argc==1){fprintf(stderr, "Usage:\n stdin | ./vcfTool COMMAND -N <sample size> [-sep <char>] [-ase] [-eigen]\n");return -1;}
	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-v")==0){verbose=atoi(argv[i+1]);}}

	for(i=0; i<argc; i++){if(strcmp(argv[i],"-log10")==0){Log10();}}
		
        long N=0;
        for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-N")==0){N=atoi(argv[i+1]);}}
	int genType=FORMAT_GT;
        if(N==0){fprintf(stderr, "Sample size is inappropriate\n"); return -1;}

	int homo=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-homo")==0){homo=1;}}

	int command=0;
	if(strcmp(argv[1],"ANALYSIS")==0){
		command=COMMAND_ANALYSIS;
	}else if(strcmp(argv[1],"CONV")==0){
		command=COMMAND_CONV;
	}else if(strcmp(argv[1],"STATS")==0){
		command=COMMAND_STATS;
	}else if(strcmp(argv[1],"FILTER")==0){
		command=COMMAND_FILTER;
	}else{
		fprintf(stderr, "No command specified!\n");
		return -1;
	}

	if(verbose>0)fprintf(stderr, "command=%d\n", command);
	int printFlag = -1;
	if(command==COMMAND_CONV){
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-dip")==0){printFlag=COMMAND_CONV_DIP;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-hapmap")==0){printFlag=COMMAND_CONV_HAPMAP;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-dose")==0){printFlag=COMMAND_CONV_DOSE;genType=FORMAT_GT;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-ase")==0){printFlag=COMMAND_CONV_AS;genType=FORMAT_GL;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-dhetratio")==0){printFlag=COMMAND_CONV_DHETR;genType=FORMAT_GL;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-gl")==0){printFlag=COMMAND_CONV_GL;genType=FORMAT_GL;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-sun")==0){printFlag=COMMAND_CONV_SUN;genType=FORMAT_GT;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-ap")==0){printFlag=COMMAND_CONV_AP;genType=FORMAT_GT;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-impute")==0){printFlag=COMMAND_CONV_IMPUTE;genType=FORMAT_GT;}}
	}else if(command==COMMAND_ANALYSIS){
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-ase")==0){printFlag=COMMAND_ANALYSIS_AS;genType=FORMAT_GL;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-gencall")==0){printFlag=COMMAND_ANALYSIS_GENCALL;genType=FORMAT_GT;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-bf")==0){printFlag=COMMAND_ANALYSIS_BF;genType=FORMAT_GT;}}
	}else if(command==COMMAND_STATS){
		printFlag=COMMAND_STATS_GENERAL;
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-rsq")==0){printFlag=COMMAND_STATS_RSQ;genType=FORMAT_GT;}}
		for(i=0; i<argc; i++){if(strcmp(argv[i],"-gencount")==0){printFlag=COMMAND_STATS_GENCOUNT;genType=FORMAT_GT;}}
	}
	int forced=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-forced")==0){forced=1;}}

	int exonOverlap=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-exonOverlap")==0){exonOverlap=1;}}
		
	double* w;
	w=(double*)calloc(N, sizeof(double));
	for(i=0; i<N; i++){w[i]=1.0;}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-w")==0){
		FILE* fw=NULL;
		fw=fopen(argv[i+1],"rb");
		if(fw==NULL){fprintf(stderr, "No proper weights\n");return -1;}
		fread(w, N, sizeof(double), fw);
		fclose(fw);
	}}
	if(verbose>10){printDouble(w,N);}
	
	
	
	// Exonic region
	char* sa=NULL;
	char* sb=NULL;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-a")==0){sa=argv[i+1];}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-b")==0){sb=argv[i+1];}}
	int nexon=0;
	long* starts;
	long* ends;
	if(sa!=NULL){
		nexon=countFields(sa);
	}
	if(nexon>0){
		starts = (long*)calloc(nexon, sizeof(long));
		ends   = (long*)calloc(nexon, sizeof(long));
		splitCSV(sa, nexon, starts);
		splitCSV(sb, nexon, ends);
		if(verbose>0)printLong(starts, nexon);
		if(verbose>0)printLong(ends, nexon);
	}
	
	// initialization for parseVCF
        //if(init()==0){fprintf(stderr, "momory allocation failuer\n"); return 0;}
	
	double MAF=-1.0;
	double IA=-1.0;
	double RSQ=-1.0;
	double HWE=1.0e16;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-MAF")==0){MAF=(double)atof(argv[i+1]);}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-IA")==0){  IA=(double)atof(argv[i+1]);}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-RSQ")==0){  RSQ=(double)atof(argv[i+1]);}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-HWE")==0){HWE=(double)atof(argv[i+1]);}}
	
	char* sep=NULL;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-sep")==0){sep=argv[i+1];}}
	if(sep==NULL){sep=(char*)calloc(2,sizeof(char));sep[0]='\t';}
	VCF_info vinfo;
	vinfo.RSQ=0.0;
        int* dip;
        double* gen;
        double* gl;
        double* ap;
	long* ase;
	double* rd;
	long pos;
	char* chr;
	char* rs;
	char* al;
        dip=(int*)calloc(N*2, sizeof(int));
        gen=(double*)calloc(N, sizeof(double));
        gl =(double*)calloc(N*3, sizeof(double));
        ap =(double*)calloc(N*2, sizeof(double));
	ase=(long*)calloc(N*2, sizeof(long));
	rd =(double*)calloc(N, sizeof(double));
	al = (char*)calloc(20, sizeof(char));
	chr=(char*)calloc(100, sizeof(char));
	rs=(char*)calloc(1000, sizeof(char));
	
	long* aseAll=NULL;
	if(command==COMMAND_ANALYSIS){
		aseAll=(long*)calloc(N*2, sizeof(long));
		clearLong(aseAll, N*2);
	}else if(command==COMMAND_STATS && printFlag==COMMAND_STATS_GENERAL){
		printf("Chrom\tPos\tRsID\tA1\tA2\tAF\tHWE\tIA\tCR\n");
	}else if(command==COMMAND_STATS && printFlag==COMMAND_STATS_RSQ){
		FILE* addfile;
		for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-rsq")==0){addfile = fopen(argv[i+1], "r");}}
		m=M=0;
		while(parseLine0(chr, &pos, rs, al, &vinfo, dip, gen, gl, ap, ase, rd, N, genType, addfile)>0){
			M++;
		}
		fclose(addfile);
		Dip = (int*)calloc(M*N*2, sizeof(int));
		rssAdd = (char**)calloc(M, sizeof(char*));
		for(m=0; m<M; m++){rssAdd[m] = (char*)calloc(100,sizeof(char)); }
		m=0;
		for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-rsq")==0){addfile = fopen(argv[i+1], "r");}}
		while(parseLine0(chr, &pos, rssAdd[m], al, &vinfo, Dip+m*N*2, gen, gl, ap, ase, rd, N, genType, addfile)>0){
                    double af = getAF(gen, w, N); //((double)tot)/((double)N)/2.0;
                    double rsq = vinfo.RSQ;
                    double ia = getIA(gl, w, N);
                    double hwe = getHWE(Dip+m*N*2, w, N);
                    if(isnan(af)>0){af=-1.0;}
                    if(isnan(ia)>0){ia=-.5;}
                    if(isnan(hwe)>0){hwe=-1.0;}
                    //if(isnan(hwe)>0){hwe=-1;}
                    if((af>MAF && af<(1.0-MAF) && isExon(pos, starts, ends, nexon)>0 && ia>IA && hwe<HWE && rsq > RSQ) || forced>0){
			m++;
		    }
                }
                M=m;
		fprintf(stderr, "Num of SNPs in the additional file=%d\n", M);
	}
	
	int numOfLine=0;
	double af;
	double ia;
	double rsq;
	double hwe;
	char sep1;
        while(parseLine(chr, &pos, rs, al, &vinfo, dip, gen, gl, ap, ase, rd, N, genType)>0){
		numOfLine++;
if(verbose>10){printDouble(gen,N);}
if(verbose>0){fprintf(stderr, "all=%c %c\n", al[0], al[1]);}
		if(command==COMMAND_STATS && exonOverlap>0){// whether snp in exon or not
			if(isExon(pos, starts, ends, nexon)>0){
				printf("%s\t1\n", rs); 
			}else{
				printf("%s\t0\n", rs);
			}
		}else{
		af = getAF(gen, w, N); //((double)tot)/((double)N)/2.0;
		rsq = vinfo.RSQ;
		ia = getIA(gl, w, N);
		hwe = getHWE(dip, w, N);
		if(isnan(af)>0){af=-1.0;}
		if(isnan(ia)>0){ia=-.5;}
		if(isnan(hwe)>0){hwe=-1.0;}
		//if(isnan(hwe)>0){hwe=-1;}
if(verbose>10){fprintf(stderr, "%s %lf %lf %lf\n", rs, af, ia, hwe);}
//fprintf(stderr, "%d %ld %ld", nexon, starts[nexon-1], ends[nexon-1]);
		if((af>MAF && af<(1.0-MAF) && isExon(pos, starts, ends, nexon)>0 && ia>IA && hwe<HWE && rsq > RSQ) || forced>0){
			if(command==COMMAND_CONV){
				if(printFlag==COMMAND_CONV_GL){
					printf("SNP%d\t%s\t%ld\t%c\t%c\t", numOfLine, rs, pos, al[0], al[1]);
				}else{
					printf("%s_%ld_%c_%c\t", chr, pos, al[0], al[1]);
				}
int aseflag=0;
int asetot=0;
				for(i=0; i<N; i++){
					sep1=sep[0];
					if(i==N-1){sep1='\n';}
					if(printFlag==COMMAND_CONV_GEN){
						printf("%1.0lf%c", gen[i],sep1);
					}else if(printFlag==COMMAND_CONV_DIP){
						printf("%d %d%c", dip[i*2],dip[i*2+1],sep1);
					}else if(printFlag==COMMAND_CONV_GL){
						printf("%lf %lf %lf%c", gl[i*3],gl[i*3+1],gl[i*3+2],sep1);
					}else if(printFlag==COMMAND_CONV_AS){
						printf("%ld %ld%c", ase[i*2],ase[i*2+1],sep1);
						//printf("%ld,%ld:%lf%c", ase[i*2],ase[i*2+1],rd[i],sep1);
						/*double ai=(double)ase[i*2]/((double)ase[i*2]+(double)ase[i*2+1]);
						if(isnan(ai)==0 && dip[i*2]+dip[i*2+1]==1 && ase[i*2]+ase[i*2+1]>5 && ai<0.99 && ai>0.01){
							aseflag++; asetot+=ase[i*2]+ase[i*2+1];
						}else{ai=-1.0;}
						if(i<N-1){
							printf("%lf%c", ai,sep1);
						}else{
							printf("%lf\t%d\t%d%c", ai, aseflag,asetot,sep1);
						}*/
						//printf("%lf%c", (double)ase[i*2]/((double)ase[i*2]+(double)ase[i*2+1]),sep1);
					}else if(printFlag==COMMAND_CONV_DOSE){
						if(i==N-1){printf("%1.0lf\n", gen[i]);}else{printf("%1.0lf", gen[i]);}
					}else if(printFlag==COMMAND_CONV_SUN){
						printf("%d%c", dip[i*2]*3+dip[i*2+1], sep1);
					}else if(printFlag==COMMAND_CONV_AP){
                                                printf("%lf\t%lf%c", ap[i*2], ap[i*2+1], sep1);
                                        }else if(printFlag==COMMAND_CONV_HAPMAP){
						printf("%c%c%c", dip[i*2]==0?al[0]:al[1], dip[i*2+1]==0?al[0]:al[1], sep1);
					}else if(printFlag==COMMAND_CONV_IMPUTE){
						if(dip[i*2]+dip[i*2+1]==0){
							printf("1 0 0%c", sep1);
						}else if(dip[i*2]+dip[i*2+1]==1){
							printf("0 1 0%c", sep1);
						}else if(dip[i*2]+dip[i*2+1]==2){
							printf("0 0 1%c", sep1);
						}else{
							printf("0 0 0%c", sep1);
						}
					}else if(printFlag==COMMAND_CONV_DHETR){
						if(ase[i*2]+ase[i*2+1]>5){
							if(dip[i*2]==1 && dip[i*2+1]==0){
								printf("%lf%c", (double)ase[i*2+1]/((double)ase[i*2]+(double)ase[i*2+1]) - gen[i]/2.0,sep1);
							}else{
								printf("%lf%c", (double)ase[i*2]/((double)ase[i*2]+(double)ase[i*2+1]) - gen[i]/2.0,sep1);
							}
						}else{
							printf("nan%c", sep1);
						}
					}
				}
			}else if(command==COMMAND_ANALYSIS){
				if(printFlag==COMMAND_ANALYSIS_AS){
					for(i=0; i<N; i++){
						int aco = dip[i*2]+dip[i*2+1];
						if(aco==1 || (homo==1 && (aco==0 || aco==2))){
							if(dip[i*2]==1){
								aseAll[i*2] += ase[i*2+1]; aseAll[i*2+1] += ase[i*2];
							}else{
								aseAll[i*2] += ase[i*2]; aseAll[i*2+1] += ase[i*2+1];
							}
						}
					}
				}else if(printFlag==COMMAND_ANALYSIS_GENCALL){
					if(covCheck(ase,N,10)>0){
						printf("%s\t%s\t%ld\t%c\t%c\t", chr, rs, pos, al[0], al[1]);
						for(i=0; i<N; i++){
							if(i<N-1){sep1='\t';}else{sep1='\n';}
							if(ase[i*2]+ase[i*2+1]>10){
								printf("%lf%c", (double)ase[i*2+1]/(double)(ase[i*2]+ase[i*2+1]), sep1);
							}else{
								printf("NA%c", sep1);
							}
						}
					}
				}else if(printFlag==COMMAND_ANALYSIS_BF){
					double thp = 0.0;
					double thn = 0.0;
					double avgbf0=0.0;
					double avgbf1=0.0;
					double top0p=0.0;
					double top0n=0.0;
					double top1p=0.0;
					double top1n=0.0;
					int imaxp=-1;
					int imaxn=-1;
					for(i=0; i<N; i++){
						avgbf0 += ap[i*2]*w[i]/(double)N;
						avgbf1 += ap[i*2+1]*w[i]/(double)N;
						if((ap[i*2]-ap[i*2+1])>thp){
						//if(fabs(ap[i*2]+ap[i*2+1])>th){
							thp = (ap[i*2]-ap[i*2+1]);//ap[i*2]>ap[i*2+1] ? ap[i*2] : ap[i*2+1];
							imaxp=i+1;
							top0p=ap[i*2];
							top1p=ap[i*2+1];
						}
						if((ap[i*2]-ap[i*2+1])<thn){
							thn = (ap[i*2]-ap[i*2+1]);
							imaxn=i+1;
							top0n=ap[i*2];
							top1n=ap[i*2+1];
						}
					}
					printf("%s\t%ld\t%s\t%c\t%c\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\n", chr, pos, rs, al[0], al[1], thp, imaxp, top0p, top1p, avgbf0, avgbf1, thn, imaxn, top0n, top1n);
				}
			}else if(command==COMMAND_STATS){
if(verbose>100){printDouble(gl,N*3);}
				if(printFlag==COMMAND_STATS_GENERAL){
					printf("%s\t%ld\t%s\t%c\t%c\t%lf\t%lf\t%lf\t%lf\n", chr, pos, rs, al[0], al[1], af, hwe, ia, getCR(gen,w,N));
				}else if(printFlag==COMMAND_STATS_RSQ){
					printf("%s\t%ld\t%s\t%c\t%c\t%lf\t%lf\t%lf\t%lf", chr, pos, rs, al[0], al[1], af, hwe, rsq, getCR(gen,w,N));
					for(m=0; m<M; m++){
						printf("\t%lf\t%s", getRsq(dip,Dip+N*m*2, N), rssAdd[m]);
						//printf("\t%lf", getRsq(dip,Dip+N*m*2, N));
					}
					printf("\n");
				}else if(printFlag==COMMAND_STATS_GENCOUNT){
					int gs[3]; gs[0]=gs[1]=gs[2]=0;
					for(i=0; i<N; i++){
						gs[dip[i*2]+dip[i*2+1]]++;
					}
					printf("%s\t%ld\t%s\t%c\t%c\t%d\t%d\t%d\n", chr, pos, rs, al[0], al[1], gs[0], gs[1], gs[2]);
				}
			}
		}}
	}
	if(command==COMMAND_ANALYSIS){
		if(printFlag==COMMAND_ANALYSIS_AS){
			//printf("%s\t%ld\t%c\t%c\t", rs, pos, al[0], al[1]);
			for(i=0; i<N-1; i++){
				printf("%ld%s%ld%s", aseAll[i*2], sep, aseAll[i*2+1], sep);
			}
			printf("%ld%s%ld\n", aseAll[i*2], sep, aseAll[i*2+1]);
		}
	}
	//print(ap,N);
        return 1;
}


