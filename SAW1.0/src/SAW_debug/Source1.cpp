#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;
int main()
{
 int i,
  j,
  k,
  y,
  m, //一维坐标
  n, //二维坐标
  l, //三维坐标
  x;

 cout << "input value for m,n,l,x:";
 cin>>m>>n>>l>>x;
 vector<vector<vector<int> > > vecInt(m, vector<vector<int> >(n, vector<int>(l)));  
 for (i = 0; i < m; i++)
  for (j = 0; j < n; j++)
   for(k = 0; k < l; k++)
	 for(y = 0;y < x; y++)
        vecInt[i][j][k][y] = i+j+k+y; 
   
 for (i = 0; i < m; i++)
 {
  for (j = 0; j < n; j++)
  {
   for(k = 0; k<l; k++)
   {
	    for(y = 0; y<x; y++)
		{
    cout<<setw(5)<<vecInt[i][j][k][y]<<":"<<setw(9)<<&vecInt[i][j][k][y];
   cout<<endl;
  }
  cout<<endl;
  }
  }
 }

 return 0;
}