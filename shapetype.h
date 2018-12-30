/* SPecification of shape*/
#ifndef SHAPETYPE_H
#define SHAPETYPE_H

#include <iostream>
#include <cmath>
using namespace std;

class ShapeType
{
public:
        ShapeType(double = 6.0, double = 5.0);

        void SetZCoordinate(double dz, int row);

        double GetXCoordinate(double theta);
        double GetYCoordinate(double theta);
        double GetZCoordinate();
        double GetHeight();
        double GetRadius();

private:
    double height, radius , Z;

    double RofTheta(double theta)
    //change R(theta) (R as a function of theta) for different
    //cross sectional geometries;
    {
         //radius+= (radius)*sin(theta)/(1-sin(theta));
        int r = 4*cos((3 * theta));
        //int r = radius;
        //double n = 10.0;

        //int r = radius / pow(pow(cos(theta),n) + pow(sin(theta),n) , 1/n);

        return r;

    }

};
#endif
//**********************Definitions****************************
ShapeType::ShapeType(double h, double r)
{
    height = h;
    radius = r;

}

void ShapeType::SetZCoordinate(double dz, int row)
{
    Z = dz*row;
}

double ShapeType::GetXCoordinate(double theta)
{


    return (RofTheta(theta) * cos(theta));
}

double ShapeType::GetYCoordinate(double theta)
{

    return (RofTheta(theta) * sin(theta)) ;
}
double ShapeType::GetZCoordinate()
{
    return Z;
}
double ShapeType::GetHeight()
{
    return height;
}
double ShapeType::GetRadius()
{
    return radius;
}
