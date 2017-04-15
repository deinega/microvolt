# include "loggerio.h"

void main() {
   char ch[3][3]={{'a', 'a', 'a'}, {'b', 'c', 'd'}, {'e', 'f', 'g'}};
   char ch2[4]={' ', ' ', ' ', ' '} ;
   
   SeqRecord S;
   S.SetSeqLen(3);
   S.OpenRecord("record.dat", "wr");
//   S.OpenRecord();
   S.AppendData(ch[0],3);
   S.AppendData(ch[0],3);
   S.NextSlice();
   S.AppendData(ch[1],3);
   S.NextSlice();
   S.GetData(ch2,2,1,4);
   S.CloseRecord();

   S.OpenRecord("record.dat", "ar");
   S.GetData(ch2,200,1,4);
   
   int i;
   for(i=0;i<100;i++){
     S.AppendData(ch[2],1);
     S.AppendData(ch[2]+1,2);
     S.NextSlice();
   }
   S.CloseRecord();

   S.SetSeqLen(100);
   S.OpenRecord("record.dat", "ar");
   

}