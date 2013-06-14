#include <fstream>
#include <vector>
#include <math.h>
#include "gauss.h"

int main()
{
  std::ifstream fin("input.txt");
	std::ofstream fout("output.txt");

	std::vector<std::vector<double> > a;
	std::vector<double> b;
	int m, n;

	fin >> m >> n;
	a.resize(m);
	b.resize(m);
	for (int i = 0; i < m; ++i)
		a[i].resize(n);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j)
			fin >> a[i][j];
		fin >> b[i];
	}

	std::vector<std::vector<double> > x;
	x = Gauss(a, b);
	if (x.empty())
		fout << "Система несовместна";
	else {
		int free_vars_count = x[0].size() - 1;
		for (int i = 0; i < n; ++i) {
			fout << "x[" << i+1 << "]" << " = " << x[i][0];
			for (int j = 1; j <= free_vars_count; ++j) {
				if (fabs(x[i][j]) > 1e-12) {
					if (x[i][j] > 0)
						fout << " + " << x[i][j] << "a[" << j << "]";
					else
						fout << " - " << -x[i][j] << "a[" << j << "]";
				}
			}
			fout << std::endl;
		}
	}

	return 0;
}
