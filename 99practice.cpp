#include<iostream>
#include<fstream>

#define F "test.txt"
#define n 9

using namespace std;

int main()
{
    int i,j,input[n][n];
    ofstream outputfile(F);

    cout<<"\n九九を表示します。\n\n"<<endl;
    outputfile<<"\n九九を表示します。\n\n\n";

    for (int i=1; i<=n; i++) {
        for (int j=1; j<=n; j++) {
            input[i][j]=i*j;
            outputfile<<input[i][j]<<"\t";
            cout<<input[i][j]<<"\t";
        }
        outputfile<<"\n";
        cout<<endl;
    }
    outputfile<<"\n";
    cout<<endl;
    outputfile.close();
}