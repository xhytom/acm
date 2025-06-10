#include<bits/stdc++.h>

using namespace std;

int main()
{
	for(int i = 1 ; i <= 10000 ; i++ )
	{
		system("E:\\cpp_files\\compare\\data.exe");
		double st = clock();
		system("E:\\cpp_files\\compare\\my.exe");
		double ed = clock();
		system("E:\\cpp_files\\compare\\std.exe");
		//此方法需要freopen();
		if(system("fc E:\\cpp_files\\compare\\my.out E:\\cpp_files\\compare\\std.out"))
		{
			puts("wrong answer");
			return 0;
		}else{
			printf("Accepted, #%d time = %.0lfms", i, ed-st);
		}
	}
}