#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static long pos[1024]; static size_t n;
static int compare_long(const void *a,const void *b){long x=*(const long*)a,y=*(const long*)b;return (x>y)-(x<y);}
static int selected(long x){return bsearch(&x,pos,n,sizeof(*pos),compare_long)!=NULL;}
int main(int argc,char **argv){
 if(argc!=2)return 2; FILE *f=fopen(argv[1],"r");if(!f)return 3;char line[8192];
 if(!fgets(line,sizeof(line),f))return 4;
 while(fgets(line,sizeof(line),f)){char *a=strchr(line,'\t');if(!a||n>=1024)return 5;pos[n++]=strtol(a+1,NULL,10);}
 if(ferror(f)||!n)return 6;fclose(f);qsort(pos,n,sizeof(*pos),compare_long);
 if(!fgets(line,sizeof(line),stdin))return 7;fputs(line,stdout);
 while(fgets(line,sizeof(line),stdin)){
  if(!strchr(line,'\n')&&!feof(stdin))return 8;
  char *a=strchr(line,'\t');if(!a)return 9;
  if(!selected(strtol(a+1,NULL,10)))continue;
  for(int i=0;i<3;i++){a=strchr(a+1,'\t');if(!a)return 10;}
  if(selected(strtol(a+1,NULL,10)))fputs(line,stdout);
 }
 if(ferror(stdin)||ferror(stdout))return 11;return 0;
}
