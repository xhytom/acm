#define PI 3.1415926
struct point//点
{
    double x,y;
};
struct circle//圆
{
    point center;
    double r;
};
double dist(point a,point b)//求圆心距
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}


double area(circle a, circle b)
{
    if((dist(a.center,b.center)+min(a.r,b.r))<=max(a.r,b.r))//两圆内交或包含
        {
            if(a.r<b.r) return PI*a.r*a.r;
            else return PI*b.r*b.r;
        }
    else if(dist(a.center,b.center)>=(a.r+b.r))//两圆外交或相离
        return 0.0;
    else
    {
        double length=dist(a.center,b.center);
        double d1=2*acos((a.r*a.r+length*length-b.r*b.r)/(2*a.r*length)); //求两扇形圆心角
        double d2=2*acos((b.r*b.r+length*length-a.r*a.r)/(2*b.r*length));
        double area1=a.r*a.r*d1/2-a.r*a.r*sin(d1)/2;
        //根据圆心角求扇形面积，减去三角形面积，得到相交部分面积
        //扇形面积：S=PI*r*r*θ/(2*PI)  三角形面积：S=1/2*a*c*sin(B)
        double area2=b.r*b.r*d2/2-b.r*b.r*sin(d2)/2;
        double area=area1+area2;
        return area;
    }
}