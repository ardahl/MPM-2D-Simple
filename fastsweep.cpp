//Testing interpolation
//g++ $(pkg-config --cflags eigen3) -std=c++11 fastsweep.cpp -o fastsweep
#include <Eigen/Dense>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <tuple>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace Eigen;

inline void fastSweep(double *field, std::vector<char> valid, double h, int res[2], int iters, double eps) {
    double err = 100;
    int it = 0;
    while(it < iters && err > eps) {
        int lowi, highi, inci;
        int lowj, highj, incj;
        switch(it % 4) {
            case 0:
                lowi = 0;
                lowj = 0;
                highi = res[0];
                highj = res[1];
                inci = 1;
                incj = 1;
                break;
            case 1:
                lowi = res[0]-1;
                lowj = 0;
                highi = -1;
                highj = res[1];
                inci = -1;
                incj = 1;
                break;
            case 2:
                lowi = res[0]-1;
                lowj = res[1]-1;
                highi = -1;
                highj = -1;
                inci = -1;
                incj = -1;
                break;
            case 3:
                lowi = 0;
                lowj = res[1]-1;
                highi = res[0];
                highj = -1;
                inci = 1;
                incj = -1;
                break;
        }
        err = 0;
        for(int i = lowi; i != highi; i += inci) {
            for(int j = lowj; j != highj; j += incj) {
                int index = i*res[1]+j;
                if(valid[index]) {
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);

                //[(x-a)+]^2 + [(x-b)+]^2 = f^2 * h^2
                //a = u_xmin
                //b = u_ymin
                double a, b;
                if(i == 0) {
                    a = field[xp1];
                }
                else if(i == res[0]-1) {
                    a = field[xm1];
                }
                else {
                    a = std::min(field[xp1], field[xm1]);
                }
                if(j == 0) {
                    b = field[yp1];
                }
                else if(j == res[0]-1) {
                    b = field[ym1];
                }
                else {
                    b = std::min(field[yp1], field[ym1]);
                }

                double f = 1;
                double fh = f*h;            //Seeing what happens with f=1
                //Eq 2.4
                double ab = a-b;
                double xbar;
                if(std::abs(ab) >= fh) {
                    xbar = std::min(a, b) + fh;
                }
                else {
                    xbar = (a+b+std::sqrt(2*f*f*h*h-ab*ab)) / 2;
                }

                double prevVal = field[index];
                field[index] = std::min(field[index], xbar);

                //Keep the max change
                double currDif = prevVal - field[index];
                if(currDif > err) {
                    err = currDif;
                }
            }
        }

        // printf("\nSweep %d\n", it+1);
        // for(int j = res[1]-1; j >= 0; j--) {
        //     for(int i = 0; i < res[0]; i++) {
        //         printf("%09.5f ", field[i*res[1]+j]);
        //     }
        //     printf("\n");
        // }
        it++;
    }
}

inline void fastSweep(Vector2d *field, std::vector<char> valid, double h, int res[2], int iters, double eps) {
    for(int dim = 0; dim < 2; dim++) {
        double err = 100;
        int it = 0;
        while(it < iters && err > eps) {
            int lowi, highi, inci;
            int lowj, highj, incj;
            switch(it % 4) {
                case 0:
                    lowi = 0;
                    lowj = 0;
                    highi = res[0];
                    highj = res[1];
                    inci = 1;
                    incj = 1;
                    break;
                case 1:
                    lowi = res[0]-1;
                    lowj = 0;
                    highi = -1;
                    highj = res[1];
                    inci = -1;
                    incj = 1;
                    break;
                case 2:
                    lowi = res[0]-1;
                    lowj = res[1]-1;
                    highi = -1;
                    highj = -1;
                    inci = -1;
                    incj = -1;
                    break;
                case 3:
                    lowi = 0;
                    lowj = res[1]-1;
                    highi = res[0];
                    highj = -1;
                    inci = 1;
                    incj = -1;
                    break;
            }
            err = 0;
            for(int i = lowi; i != highi; i += inci) {
                for(int j = lowj; j != highj; j += incj) {
                    int index = i*res[1]+j;
                    if(valid[index]) {
                        continue;
                    }
                    int xp1 = (i+1)*res[1]+j;
                    int xm1 = (i-1)*res[1]+j;
                    int yp1 = i*res[1]+(j+1);
                    int ym1 = i*res[1]+(j-1);

                    //[(x-a)+]^2 + [(x-b)+]^2 = f^2 * h^2
                    //a = u_xmin
                    //b = u_ymin
                    double a, b;
                    if(i == 0) {
                        a = field[xp1](dim);
                    }
                    else if(i == res[0]-1) {
                        a = field[xm1](dim);
                    }
                    else {
                        a = std::min(field[xp1](dim), field[xm1](dim));
                    }
                    if(j == 0) {
                        b = field[yp1](dim);
                    }
                    else if(j == res[0]-1) {
                        b = field[ym1](dim);
                    }
                    else {
                        b = std::min(field[yp1](dim), field[ym1](dim));
                    }

                    double f = 1;
                    double fh = f*h;            //Seeing what happens with f=1
                    //Eq 2.4
                    double ab = a-b;
                    double xbar;
                    if(std::abs(ab) >= fh) {
                        xbar = std::min(a, b) + fh;
                    }
                    else {
                        xbar = (a+b+std::sqrt(2*f*f*h*h-ab*ab)) / 2;
                    }

                    double prevVal = field[index](dim);
                    field[index](dim) = std::min(field[index](dim), xbar);

                    //Keep the max change
                    double currDif = prevVal - field[index](dim);
                    if(currDif > err) {
                        err = currDif;
                    }
                }
            }

            // printf("\nSweep %d, Dim %d\n", it+1, dim);
            // for(int j = res[1]-1; j >= 0; j--) {
            //     for(int i = 0; i < res[0]; i++) {
            //         printf("(%09.5f,%09.5f) ", field[i*res[1]+j](0), field[i*res[1]+j](1));
            //     }
            //     printf("\n");
            // }
            it++;
        }
    }
}

inline void velExtrapolateFM(Vector2d *vel, double *field, std::vector<char> valid, double h, int res[2]) {
    //Sort phi by smallest to largest
    std::vector<std::tuple<double,int,int>> sortedphi;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            if(!valid[i*res[1]+j]) {
                //Could switch but then sort requires custom function
                sortedphi.push_back(std::make_tuple(field[i*res[1]+j], i, j));
            }
        }
    }
    std::sort(sortedphi.begin(), sortedphi.end());

    //Go through array, calculating new velocity for each element
    std::vector<std::tuple<double,int,int>>::iterator it;
    bool d1, d2;
    double tphi, phi1, phi2;
    Vector2d x1, x2;
    for(it = sortedphi.begin(); it != sortedphi.end(); it++) {
        double phi = std::get<0>(*it);
        d1 = d2 = false;
        int x = std::get<1>(*it);
        int y = std::get<2>(*it);
        int index = x*res[1]+y;
        int xp1 = (x+1)*res[1]+y;
        int xm1 = (x-1)*res[1]+y;
        int yp1 = x*res[1]+(y+1);
        int ym1 = x*res[1]+(y-1);

        //Look in x direction
        if(x != 0 && field[xm1] < phi) {
            phi1 = field[xm1];
            d1 = true;
            x1 = vel[xm1];
        }
        if(x != res[0]-1 && field[xp1] < phi) {
            if((d1 && field[xp1] < phi1) || !d1) {
                phi1 = field[xp1];
                d1 = true;
                x1 = vel[xp1];
            }
        }

        //look in y direction
        if(y != 0 && field[ym1] < phi) {
            phi2 = field[ym1];
            d2 = true;
            x2 = vel[ym1];
        }
        if(y != res[1]-1 && field[yp1] < phi) {
            if((d2 && field[yp1] < phi2) || !d2) {
                phi2 = field[yp1];
                d2 = true;
                x2 = vel[yp1];
            }
        }

        //Compute average value
        if(d1 && d2) {
            vel[index] = (x1*(phi-phi1)+x2*(phi-phi2)) / (2*phi-phi1-phi2);
        }
        else if(d1) {
            vel[index] = x1;
        }
        else if(d2) {
            vel[index] = x2;
        }
        else {
            vel[index] = Vector2d(0,0);
        }
    }
}

inline void velExtrapolateFS(Vector2d *vel, double *field, std::vector<char> valid, int res[2], int iters, double eps) {
    double err = 100;
    int it = 0;
    while(it < iters && err > eps) {
        int lowi, highi, inci;
        int lowj, highj, incj;
        switch(it % 4) {
            case 0:
                lowi = 0;
                lowj = 0;
                highi = res[0];
                highj = res[1];
                inci = 1;
                incj = 1;
                break;
            case 1:
                lowi = res[0]-1;
                lowj = 0;
                highi = -1;
                highj = res[1];
                inci = -1;
                incj = 1;
                break;
            case 2:
                lowi = res[0]-1;
                lowj = res[1]-1;
                highi = -1;
                highj = -1;
                inci = -1;
                incj = -1;
                break;
            case 3:
                lowi = 0;
                lowj = res[1]-1;
                highi = res[0];
                highj = -1;
                inci = 1;
                incj = -1;
                break;
        }
        err = 0;
        for(int i = lowi; i != highi; i += inci) {
            for(int j = lowj; j != highj; j += incj) {
                int index = i*res[1]+j;
                if(valid[index]) {
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);

                //(Fxmin*(phix-phixmin) + Fymin*(phix-phiymin)) / ((phix-phixmin) + (phix-phiymin))
                double a, b;
                Vector2d vi, vj;
                if(i == 0) {
                    a = field[xp1];
                    vi = vel[xp1];
                }
                else if(i == res[0]-1) {
                    a = field[xm1];
                    vi = vel[xm1];
                }
                else {
                    if(field[xp1] < field[xm1]) {
                        a = field[xp1];
                        vi = vel[xp1];
                    }
                    else {
                        a = field[xm1];
                        vi = vel[xm1];
                    }
                }
                if(j == 0) {
                    b = field[yp1];
                    vj = vel[yp1];
                }
                else if(j == res[0]-1) {
                    b = field[ym1];
                    vj = vel[ym1];
                }
                else {
                    if(field[yp1] < field[ym1]) {
                        b = field[yp1];
                        vj = vel[yp1];
                    }
                    else {
                        b = field[ym1];
                        vj = vel[ym1];
                    }
                }

                //If neither values are less than x_ij then the contribution is 0
                double phixi = 0, phixj = 0;
                if(field[index] > a) {
                    phixi = field[index]-a;
                }
                if(field[index] > b) {
                    phixj = field[index]-b;
                }
                //Solving eq. at bottom of pg 13 in Adalsteinsson and Sethian 1999
                Vector2d vtmp = (vi*phixi + vj*phixj) / (phixi + phixj);
                double dnorm = (vel[index]-vtmp).norm();
                if(dnorm > err) {
                    err = dnorm;
                }
                vel[index] = vtmp;
            }
        }
        // printf("\nSweep %d\n", it+1);
        // for(int j = res[1]-1; j >= 0; j--) {
        //     for(int i = 0; i < res[0]; i++) {
        //         printf("(%09.5f,%09.5f) ", vel[i*res[1]+j](0), vel[i*res[1]+j](1));
        //     }
        //     printf("\n");
        // }
        it++;
    }
}

inline void velExtrapolateAve(Vector2d* vel, std::vector<char> valid, int res[2], double iters, double eps) {
    std::vector<char> oldValid;
    oldValid.resize(res[0]*res[1]);
    Vector2d *tmpVel = new Vector2d[res[0]*res[1]];
    //Perform this a couple of times
    int it = 0;
    double err = 100;
    while(it < iters && err > eps) {
    /// for(int it = 0; it < iter; it++) {
        oldValid = valid;
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);

                Vector2d sum(0,0);
                int count = 0;

                tmpVel[ind] = vel[ind];
                //Only check around cells that are not already valid
                if(!oldValid[ind]) {
                    if(i+1 < res[0] && oldValid[xp1]) {
                        sum += vel[xp1];
                        count++;
                    }
                    if(i-1 >= 0 && oldValid[xm1]) {
                        sum += vel[xm1];
                        count++;
                    }
                    if(j+1 < res[1] && oldValid[yp1]) {
                        sum += vel[yp1];
                        count++;
                    }
                    if(j-1 >= 0 && oldValid[ym1]) {
                        sum += vel[ym1];
                        count++;
                    }
                    if(i+1 < res[0] && j+1 < res[1] && oldValid[(i+1)*res[1]+(j+1)]) {
                        sum += vel[(i+1)*res[1]+(j+1)];
                        count++;
                    }
                    if(i+1 < res[0] && j-1 >=0 && oldValid[(i+1)*res[1]+(j-1)]) {
                        sum += vel[(i+1)*res[1]+(j-1)];
                        count++;
                    }
                    if(i-1 >= 0 && j+1 < res[1] && oldValid[(i-1)*res[1]+(j+1)]) {
                        sum += vel[(i-1)*res[1]+(j+1)];
                        count++;
                    }
                    if(i-1 >= 0 && j-1 >= 0 && oldValid[(i-1)*res[1]+(j-1)]) {
                        sum += vel[(i-1)*res[1]+(j-1)];
                        count++;
                    }

                    //If neighboring cells were valid, average the values
                    if(count > 0) {
                        tmpVel[ind] = sum / count;
                        valid[ind] = 1;
                    }
                }
            }
        }
        err = 0;
        //Replace grid with temp grid
        for(int i = 0; i < res[0]; i++) {
            for(int j = 0; j < res[1]; j++) {
                int ind = i*res[1]+j;
                double dnorm = (vel[ind]-tmpVel[ind]).norm();
                if(dnorm > err) {
                    err = dnorm;
                }
                vel[i*res[1]+j] = tmpVel[i*res[1]+j];
            }
        }

        // printf("\nSweep %d\n", it+1);
        // for(int j = res[1]-1; j >= 0; j--) {
        //     for(int i = 0; i < res[0]; i++) {
        //         printf("(%09.5f,%09.5f) ", vel[i*res[1]+j](0), vel[i*res[1]+j](1));
        //     }
        //     printf("\n");
        // }
        it++;
    }
}

inline void matExtrapolateFM(Vector2d* mats, double *phi, std::vector<char>& valid, double h, int res[2]) {
	std::vector<std::tuple<double,int,int>> sortedphi;
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            if(!valid[i*res[1]+j]) {
                //Could switch but then sort requires custom function
                sortedphi.push_back(std::make_tuple(phi[i*res[1]+j], i, j));
            }
        }
    }
    std::sort(sortedphi.begin(), sortedphi.end());

    //Go through array, calculating new velocity for each element
    std::vector<std::tuple<double,int,int>>::iterator it;
	double xh = 0, yh = 0;
	Vector2d xm, ym;
    for(it = sortedphi.begin(); it != sortedphi.end(); it++) {
        int x = std::get<1>(*it);
        int y = std::get<2>(*it);
        int index = x*res[1]+y;
        int xp1 = (x+1)*res[1]+y;
        int xm1 = (x-1)*res[1]+y;
        int yp1 = x*res[1]+(y+1);
        int ym1 = x*res[1]+(y-1);

        //Look in x direction
        if(x < res[0]-1 && valid[xp1]) {	//right valid, -h
			xh = -h;
			xm = mats[xp1];
		}
		else if(x > 0 && valid[xm1]) {		//left valid, +h
			xh = h;
			xm = mats[xm1];
		}
		else {								//shouldn't happen, 0
			xh = 0;
			xm = Vector2d::Zero();
		}

        //look in y direction
        if(y < res[1]-1 && valid[yp1]) {	//right valid, -h
			yh = -h;
			ym = mats[yp1];
		}
		else if(y > 0 && valid[ym1]) {		//left valid, +h
			yh = h;
			ym = mats[ym1];
		}
		else {								//shouldn't happen, 0
			yh = 0;
			ym = Vector2d::Zero();
		}

		//If one direction is empty, we want to keep the value from the other direction
		if(xh == 0 && yh != 0) {
			xm = ym;
		}
		else if(yh == 0 && xh != 0) {
			ym = xm;
		}
		else if(xh == 0 && yh == 0){
			printf("Really shouldn't happen");
		}
		mats[index] = Vector2d(xm(0)+xh, ym(1)+yh);
		valid[index] = 1;
    }
}

inline void matExtrapolateFS(Vector2d *mats, double *phi, Matrix2d F, std::vector<char> valid, double h, int res[2], int iters, double eps) {
    double err = 100;
    int it = 0;
    while(it < iters && err > eps) {
        int lowi, highi, inci;
        int lowj, highj, incj;
        switch(it % 4) {
            case 0:
                lowi = 0;
                lowj = 0;
                highi = res[0];
                highj = res[1];
                inci = 1;
                incj = 1;
                break;
            case 1:
                lowi = res[0]-1;
                lowj = 0;
                highi = -1;
                highj = res[1];
                inci = -1;
                incj = 1;
                break;
            case 2:
                lowi = res[0]-1;
                lowj = res[1]-1;
                highi = -1;
                highj = -1;
                inci = -1;
                incj = -1;
                break;
            case 3:
                lowi = 0;
                lowj = res[1]-1;
                highi = res[0];
                highj = -1;
                inci = 1;
                incj = -1;
                break;
        }
        err = 0;
        for(int i = lowi; i != highi; i += inci) {
            for(int j = lowj; j != highj; j += incj) {
                int index = i*res[1]+j;
                if(valid[index]) {
                    continue;
                }
                int xp1 = (i+1)*res[1]+j;
                int xm1 = (i-1)*res[1]+j;
                int yp1 = i*res[1]+(j+1);
                int ym1 = i*res[1]+(j-1);

                //Look towards the object, direction with smallest phi
                double p = phi[index];
                double phix = p, phiy = p;
                bool positivex = true, positivey = true;
					//x direction
				if(i != res[0]-1 && phi[xp1] < phix) {
					phix = phi[xp1];
				}
				else if(i != 0 && phi[xm1] < phix) {
					phix = phi[xm1];
					positivex = false;
				}
				if(j != res[1]-1 && phi[yp1] < phiy) {
					phiy = phi[yp1];
				}
				else if(j != 0 && phi[ym1] < phiy) {
					phiy = phi[ym1];
					positivey = false;
				}
                //Solve grad(X)=F
                double x1, x2, y1, y2;
                if(positivex) {
					x1 = mats[xp1](0) - h*F(0,0);
					x2 = mats[xp1](1) - h*F(0,1);
				}
				else {
					x1 = h*F(0,0) + mats[xm1](0);
					x2 = h*F(0,1) + mats[xm1](1);
				}
				if(positivey) {
					y1 = mats[yp1](0) - h*F(1,0);
					y2 = mats[yp1](1) - h*F(1,1);
				}
				else {
					y1 = h*F(1,0) + mats[yp1](0);
					y2 = h*F(1,1) + mats[yp1](1);
				}
                //Take only the smaller one?
                //Average?
                Vector2d npos((x1+x2)/2.0, (y1+y2)/2.0);
                //Repeat until convergence
                double norm = (npos - mats[index]).norm();
				if(norm > err) {
					err = norm;
				}
                mats[index] = npos;
            }
        }
        /// printf("\nSweep %d\n", it+1);
        /// for(int j = res[1]-1; j >= 0; j--) {
            /// for(int i = 0; i < res[0]; i++) {
                /// printf("(%09.5f,%09.5f) ", vel[i*res[1]+j](0), vel[i*res[1]+j](1));
            /// }
            /// printf("\n");
        /// }
        it++;
    }
}

inline Matrix2d upwindJac(Vector2d* field, Vector2d* vel, int i, int j, std::vector<char> valid, double h, int res[2], bool boundary = true) {
    //TODO: Points inside the object don't look outside and points outside dont look inside
    //Compute D(field) with upwinding scheme
    int index = i*res[1]+j;
    int xp1 = (i+1)*res[1]+j;
    int xm1 = (i-1)*res[1]+j;
    int yp1 = i*res[1]+(j+1);
    int ym1 = i*res[1]+(j-1);
    //DX
    int di = 100, dj = 100;
    //Dv
    /// int di = 31, dj = 48;
    if(i == di && j == dj) {
		printf("(i, j): (%f, %f), index: %d\n", field[index](0), field[index](1), (int)valid[index]);
		printf("(i+1, j): (%f, %f), xp1: %d\n", field[xp1](0), field[xp1](1), (int)valid[xp1]);
		printf("(i-1, j): (%f, %f), xm1: %d\n", field[xm1](0), field[xm1](1), (int)valid[xm1]);
		printf("(i, j+1): (%f, %f), yp1: %d\n", field[yp1](0), field[yp1](1), (int)valid[yp1]);
		printf("(i, j-1): (%f, %f), ym1: %d\n", field[ym1](0), field[ym1](1), (int)valid[ym1]);
	}
    double xx = 0, xy = 0, yx = 0, yy = 0;      //Upwind scheme so look at each value individually
    if(vel[index](0) >= 0) {    //X Velocity is positive, so moving from left to right, look left
		if(i == di && j == dj) {
			printf("Look X -1\n");
		}
        if((i != 0) && (((boundary && valid[index] && valid[xm1]) || !boundary) || ((boundary && !valid[index] && !valid[xm1]) || !boundary))) {  //Looking backwards in x and valid value
			xx = (field[index](0) - field[xm1](0)) / h;
			if(i == di && j == dj) {
				printf("xx1: (%f - %f) / %f = %f\n", field[index](0), field[xm1](0), h, xx);
			}
        }
        else if((i != res[0]-1) && (((boundary && valid[index] && valid[xp1]) || !boundary) || ((boundary && !valid[index] && !valid[xp1]) || !boundary))) {  //Can't look backwards in x so do forward difference instead
			xx = (field[xp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xx2: (%f - %f) / %f = %f\n", field[xp1](0), field[index](0), h, xx);
			}
        }

        if((j != 0) && (((boundary && valid[index] && valid[ym1]) || !boundary) || ((boundary && !valid[index] && !valid[ym1]) || !boundary))) {  //Looking backwards in y and valid value
			xy = (field[index](0) - field[ym1](0)) / h;
			if(i == di && j == dj) {
				printf("xy1: (%f - %f) / %f = %f\n", field[index](0), field[ym1](0), h, xy);
			}
        }
        else if((j != res[1]-1) && (((boundary && valid[index] && valid[yp1]) || !boundary) || ((boundary && !valid[index] && !valid[yp1]) || !boundary))) {  //Can't look backwards in y so do forward difference instead
			xy = (field[yp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xy2: (%f - %f) / %f = %f\n", field[yp1](0), field[index](0), h, xy);
			}
        }
    }
    else {                          //X Velocity is negative, so forward difference
		if(i == di && j == dj) {
			printf("Look X +1\n");
		}
        if((i != res[0]-1) && (((boundary && valid[index] && valid[xp1]) || !boundary) || ((boundary && !valid[index] && !valid[xp1]) || !boundary))) {
			xx = (field[xp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xx3: (%f - %f) / %f = %f\n", field[xp1](0), field[index](0), h, xx);
			}
        }
        else if((i != 0) && (((boundary && valid[index] && valid[xm1]) || !boundary) || ((boundary && !valid[index] && !valid[xm1]) || !boundary))) {
			xx = (field[index](0) - field[xm1](0)) / h;
			if(i == di && j == dj) {
				printf("xx4: (%f - %f) / %f = %f\n", field[index](0), field[xm1](0), h, xx);
			}
        }

        if((j != res[1]-1) && (((boundary && valid[index] && valid[yp1]) || !boundary) || ((boundary && !valid[index] && !valid[yp1]) || !boundary))) {
			xy = (field[yp1](0) - field[index](0)) / h;
			if(i == di && j == dj) {
				printf("xy3: (%f - %f) / %f = %f\n", field[yp1](0), field[index](0), h, xy);
			}
        }
        else if((j != 0) && (((boundary && valid[index] && valid[ym1]) || !boundary) || ((boundary && !valid[index] && !valid[ym1]) || !boundary))) {
			xy = (field[index](0) - field[ym1](0)) / h;
			if(i == di && j == dj) {
				printf("xy4: (%f - %f) / %f = %f\n", field[index](0), field[ym1](0), h, xy);
			}
        }
    }
    if(vel[index](1) >= 0) {    //Y Velocity is positive, so backward difference
		if(i == di && j == dj) {
			printf("Look Y -1\n");
		}
        if((i != 0) && (((boundary && valid[index] && valid[xm1]) || !boundary) || ((boundary && !valid[index] && !valid[xm1]) || !boundary))) {
			yx = (field[index](1) - field[xm1](1)) / h;
			if(i == di && j == dj) {
				printf("yx1: (%f - %f) / %f = %f\n", field[index](1), field[xm1](1), h, yx);
			}
        }
        else if((i != res[0]-1) && (((boundary && valid[index] && valid[xp1]) || !boundary) || ((boundary && !valid[index] && !valid[xp1]) || !boundary))) {
			yx = (field[xp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yx2: (%f - %f) / %f = %f\n", field[xp1](1), field[index](1), h, yx);
			}
        }

        if((j != 0) && (((boundary && valid[index] && valid[ym1]) || !boundary) || ((boundary && !valid[index] && !valid[ym1]) || !boundary))) {
			yy = (field[index](1) - field[ym1](1)) / h;
			if(i == di && j == dj) {
				printf("yy1: (%f - %f) / %f = %f\n", field[index](1), field[ym1](1), h, yy);
			}
        }
        else if((j != res[1]-1) && (((boundary && valid[index] && valid[yp1]) || !boundary) || ((boundary && !valid[index] && !valid[yp1]) || !boundary))) {
			yy = (field[yp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yy2: (%f - %f) / %f = %f\n", field[yp1](1), field[index](1), h, yy);
			}
        }
    }
    else {                          //Y Velocity is negative, so forward difference
		if(i == di && j == dj) {
			printf("Look Y +1\n");
		}
        if((i != res[0]-1) && (((boundary && valid[index] && valid[xp1]) || !boundary) || ((boundary && !valid[index] && !valid[xp1]) || !boundary))) {
			yx = (field[xp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yx3: (%f - %f) / %f = %f\n", field[xp1](1), field[index](1), h, yx);
			}
        }
        else if((i != 0) && (((boundary && valid[index] && valid[xm1]) || !boundary) || ((boundary && !valid[index] && !valid[xm1]) || !boundary))) {
			yx = (field[index](1) - field[xm1](1)) / h;
			if(i == di && j == dj) {
				printf("yx4: (%f - %f) / %f = %f\n", field[index](1), field[xm1](1), h, yx);
			}
        }

        if((j != res[1]-1) && (((boundary && valid[index] && valid[yp1]) || !boundary) || ((boundary && !valid[index] && !valid[yp1]) || !boundary))) {
			yy = (field[yp1](1) - field[index](1)) / h;
			if(i == di && j == dj) {
				printf("yy3: (%f - %f) / %f = %f\n", field[yp1](1), field[index](1), h, yy);
			}
        }
        else if((j != 0) && (((boundary && valid[index] && valid[ym1]) || !boundary) || ((boundary && !valid[index] && !valid[ym1]) || !boundary))) {
			yy = (field[index](1) - field[ym1](1)) / h;
			if(i == di && j == dj) {
				printf("yy4: (%f - %f) / %f = %f\n", field[index](1), field[ym1](1), h, yy);
			}
        }
    }
    Matrix2d D;
    D << xx, xy, yx, yy;
    return D;
}

int main(int argc, char** argv) {
    int res[2] = {15, 15};
    Vector2d origin(-7.5, 0);
    double h = 1;
    Vector2d center = origin + h*Vector2d(res[0]-1, res[1]-1)/2.0;

    Vector2d *material = new Vector2d[res[0]*res[1]];
    Vector2d *field = new Vector2d[res[0]*res[1]];
    Vector2d *velExact = new Vector2d[res[0]*res[1]];
    double *distField = new double[res[0]*res[1]];
    double *distExact = new double[res[0]*res[1]];
    Matrix2d *DX = new Matrix2d[res[0]*res[1]];
    std::vector<char> valid;
    Vector2d *exact = new Vector2d[res[0]*res[1]];
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            field[i*res[1]+j] = Vector2d(100,100);
            distField[i*res[1]+j] = 100;
            valid.push_back(0);
        }
    }

    double rad = 2.5;
    for(int i = 0; i < res[0]; i++) {
		for(int j = 0; j < res[1]; j++) {
			int index = i*res[1]+j;
			Vector2d p = origin+h*Vector2d(i, j);
			Vector2d ph = p-center;
			Rotation2D<double> rot(-0.785398);
			Vector2d rotate = rot * ph;
			/// exact[index] = rotate;
			velExact[index] = 0.5*Vector2d(-ph(1), ph(0));
			exact[index] = p;
			double d = ph.norm() - rad;
			distExact[index] = d;
			material[index] = rotate+center;
			if(d < 0.0) {
                Vector2d np = rotate+center;
				field[index] = np - p;
				// field[index] = p;
				// distField[index] = d;
                distField[index] = -h/2.0;
				valid[index] = 1;
			}
		}
	}

	// std::vector<char> validcpy = valid;
	// matExtrapolateFM(material, distExact, valid, h, res);
	// /// Rotation2D<double> rot(-0.785398);
	// /// matExtrapolateFS(material, distExact, rot.toRotationMatrix(), valid, h, res, 20, 1e-8);
	// std::ofstream matout("material");
	// printf("Material Field:\n");
	// double err = 0;
    // double maxd = 0;
    // Vector2d md1, md2;
    // int mi, mj, mi2, mj2;
	// for(int j = res[1]-1; j >= 0; j--) {
	// 	for(int i = 0; i < res[0]; i++) {
	// 		int ind = i*res[1]+j;
	// 		Matrix2d dx = upwindJac(material, velExact, i, j, valid, h, res);
	// 		DX[ind] = dx;
	// 		Vector2d p = origin + h*Vector2d(i, j);
	// 		Vector2d ph = p-center;
	// 		Rotation2D<double> rot(-0.785398);
	// 		Vector2d rotate = rot * ph;
	// 		p = rotate+center;
	// 		Vector2d m = material[ind];
	// 		if(validcpy[ind]) {
	// 			printf("\033[0;31m(% 04.1f,% 04.1f)\033[0m ", m(0), m(1));
	// 		}
	// 		else {
	// 			printf("(% 04.1f,% 04.1f) ", m(0), m(1));
	// 		}
	// 		matout << p(0) << " " << p(1) << " " << "0\n";
	// 		matout << m(0) << " " << m(1) << " " << "1\n\n\n";
	// 		double diff = (p-m).norm();
	// 		if(diff > err) {
	// 			err = diff;
	// 		}
    //         if(i != res[0]-1) {
    //             Vector2d mtmp = material[(i+1)*res[1]+j];
    //             diff = (m-mtmp).norm();
    //             if(diff > maxd) {
    //                 maxd = diff;
    //                 mi = i;
    //                 mj = j;
    //                 mi2 = i+1;
    //                 mj2 = j;
    //                 md1 = m;
    //                 md2 = mtmp;
    //             }
    //         }
    //         if(i != 0) {
    //             Vector2d mtmp = material[(i-1)*res[1]+j];
    //             diff = (m-mtmp).norm();
    //             if(diff > maxd) {
    //                 maxd = diff;
    //                 mi = i;
    //                 mj = j;
    //                 mi2 = i-1;
    //                 mj2 = j;
    //                 md1 = m;
    //                 md2 = mtmp;
    //             }
    //         }
    //         if(j != res[1]-1) {
    //             Vector2d mtmp = material[i*res[1]+(j+1)];
    //             diff = (m-mtmp).norm();
    //             if(diff > maxd) {
    //                 maxd = diff;
    //                 mi = i;
    //                 mj = j;
    //                 mi2 = i;
    //                 mj2 = j+1;
    //                 md1 = m;
    //                 md2 = mtmp;
    //             }
    //         }
    //         if(j != 0) {
    //             Vector2d mtmp = material[i*res[1]+(j-1)];
    //             diff = (m-mtmp).norm();
    //             if(diff > maxd) {
    //                 maxd = diff;
    //                 mi = i;
    //                 mj = j;
    //                 mi2 = i;
    //                 mj2 = j-1;
    //                 md1 = m;
    //                 md2 = mtmp;
    //             }
    //         }
	// 	}
	// 	printf("\n");
	// }
	// printf("Max Error: %f\n", err);
	// matout.close();
    // printf("Max Dist: %f\n", maxd);
    // printf("%d, %d: (%f, %f)\n", mi, mj, md1(0), md1(1));
    // printf("%d, %d: (%f, %f)\n", mi2, mj2, md2(0), md2(1));
    //
	// std::ofstream debug("fsmat");
	// debug << "\nDX: \n";
	// for(int j = res[1]-2; j >= 0; j--) {
	// 	for(int k = 0; k < 2; k++) {
	// 		for(int i = 0; i < res[0]-1; i++) {
	// 			int ind = i*res[0]+j;
	// 			debug << "[";
	// 			if(DX[ind](k,0) < 0) {
	// 				debug << std::fixed << std::setprecision(6) << DX[ind](k,0);
	// 			}
	// 			else {
	// 				debug << std::fixed << std::setprecision(6) << DX[ind](k,0);
	// 			}
	// 			debug << ",";
	// 			if(DX[ind](k,1) < 0) {
	// 				debug << std::fixed << std::setprecision(6) << DX[ind](k,1);
	// 			}
	// 			else {
	// 				debug << std::fixed << std::setprecision(6) << DX[ind](k,1);
	// 			}
	// 			debug << "] ";
	// 		}
	// 		debug << "\n";
	// 	}
	// 	debug << "\n";
	// }
	// debug << "\n";
	// debug.close();
    // std::exit(0);

    std::vector<char> tmpvalid = valid;
    /// for(int i = 0; i < res[0]; i++) {
		/// for(int j = 0; j < res[1]; j++) {
			/// int index = i*res[1]+j;
			/// if(valid[index]) {
				/// continue;
			/// }
			/// int xp1 = (i+1)*res[1]+j;
			/// int xm1 = (i-1)*res[1]+j;
			/// int yp1 = i*res[1]+(j+1);
			/// int ym1 = i*res[1]+(j-1);
			/// if((i != res[0]-1 && valid[xp1]) || (i != 0 && valid[xm1]) ||
			   /// (j != res[1]-1 && valid[yp1]) || (j != 0 && valid[ym1])) {
				/// Vector2d p = origin+h*Vector2d(i, j);
				/// double d = (p-center).norm() - rad;
				/// distField[index] = d;
				/// tmpvalid[index] = 1;
			/// }
		/// }
    /// }

    printf("Velocity Field\n");
    for(int j = res[1]-1; j >= 0; j--) {
        for(int i = 0; i < res[0]; i++) {
            printf("(%f, %f) ", field[i*res[1]+j](0), field[i*res[1]+j](1));
        }
        printf("\n");
    }
    printf("Distance Field\n");
    for(int j = res[1]-1; j >= 0; j--) {
        for(int i = 0; i < res[0]; i++) {
            printf("%09.5f ", distField[i*res[1]+j]);
        }
        printf("\n");
    }
    /// std::exit(0);

    printf("\n\n");
    printf("----------------\n");
    printf("Distance Field\n");
    printf("----------------\n\n");
    fastSweep(distField, tmpvalid, h, res, 6, 1e-8);
    printf("\nResult\n");
    for(int j = res[1]-1; j >= 0; j--) {
        for(int i = 0; i < res[0]; i++) {
            printf("%09.5f ", distField[i*res[1]+j]);
        }
        printf("\n");
    }

    printf("\n\n");
    printf("----------------\n");
    printf("Velocity Field\n");
    printf("----------------\n\n");
    /// fastSweep(field, valid, h, res, 6, 1e-8);
    /// velExtrapolateFM(field, distField, valid, h, res);
    velExtrapolateFS(field, distField, valid, res, 10, 1e-8);
    /// velExtrapolateAve(field, valid, res, 10, 1e-8);
    printf("\nResult\n");
    for(int j = res[1]-1; j >= 0; j--) {
        for(int i = 0; i < res[0]; i++) {
            printf("(%09.5f,%09.5f) ", field[i*res[1]+j](0), field[i*res[1]+j](1));
        }
        printf("\n");
    }

    std::ofstream diff("veldiff.txt");
    std::ofstream disp("displacement.txt");
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            int index = i*res[1]+j;
            Vector2d xg = origin+h*Vector2d(i,j);
            Vector2d ph = xg-center;
            Rotation2D<double> rot(-0.785398);
            Vector2d rotate = rot * ph;
            Vector2d np = rotate+center;
            if(valid[index]) {
                disp << xg(0) << " " << xg(1) << " " << (field[index]-xg)(0) << " " <<  (field[index]-xg)(1) << " 0 255 0\n";
            }
            // else {
            //     disp << xg(0) << " " << xg(1) << " " << (field[index]-xg)(0) << " " <<  (field[index]-xg)(1) << " 255 0 0\n";
            // }
            if(field[index](0) < 100 && field[index](1) < 100) {
                Vector2d d = field[index]-np;
                diff << xg(0) << " " << xg(1) << " " << d.norm() << "\n";
            }
        }
    }

    return 0;
}
