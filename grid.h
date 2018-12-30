#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <cstdlib> //malloc
#include <fstream>

using namespace std;
//A N by M parametization of a closed surface
class Grid {

    friend void Print(ofstream&, Grid&);

public:
    Grid(int = 5, int = 5);
    ~Grid()
    {
        //deallocate memory block
        free (triangle_vertices);
    }

    int GetN();
    int GetM();
    int GetRectangle();
    int GetTriangle();
    int GetNumberofVertices();
    int GetMCoordinate(int);
    int GetNCoordinate(int);

    void SetTriangleVertices();

    int* triangle_vertices;

private:
    //helper functions
    void SecondTriangleVertices(int, int, int, int&, int&, int&);
    void AddToArray(int, int, int);
    //n(vertical, m(horizontal) grid dimensions
    int n, m, number_of_vertices, m_coordinate , n_coordinate;
    int rectangle;
    int triangle;
};
#endif
//**************** Function Definitions*******************************
Grid::Grid(int length, int width)
{
    n = length; m = width;
    rectangle = n * m;
    triangle = rectangle * 2;
    number_of_vertices = m * (n+1);
    //Dynamic 2D array, 3 vertices per triangle
    triangle_vertices = (int*) malloc( sizeof(int) * triangle * 3);
}
void Grid::SetTriangleVertices()
{
    int rectangle_node = 0;    //bottom left vertex
    int rectangle = 0;         // loops through all rectangles in grid
    int last_rectangle = m;    //last rectangle in the row
    int num_rectangles = m*n;  //last rectangle in grid
    int a ,b,c;
    //initialize to first node;
    a = 0; b = 1; c = m;

    while(rectangle < (num_rectangles))
    {
        if(a == 0)
        //at first rectangle in grid
            AddToArray(a,b,c);
        else{
            if(((last_rectangle/rectangle_node) == 1) &&
                (last_rectangle % rectangle_node) == 1)
                //If at the second to last rectangle in the row
                {
                    b = a - (m - 1);
                    c++;
                    last_rectangle += m;  //updates last rectangle for the next row
                }
                else {
                    b = a+1;
                    c++;
                }
                AddToArray(a,b,c);
            }
        a++;
        //loop control
        rectangle_node++;
        rectangle++;
        // variables represent the same concept
        //rewrite
    }
}
void Grid::SecondTriangleVertices(int i1, int j1, int k1,
                                  int& i2, int& j2, int& k2)
{
    i2 = k1;
    j2 = j1;
    k2 = j2 + m;
}

void Grid::AddToArray(int i1, int j1, int k1)
{
    static int area = 1;
    int i2, j2, k2;

    triangle_vertices[3*(area-1)]     = i1;
    triangle_vertices[3*(area-1) + 1] = j1;
    triangle_vertices[3*(area-1) + 2] = k1;

    area++;

    SecondTriangleVertices( i1,j1,k1,i2,j2,k2);

    triangle_vertices[3*(area-1)]     = i2;
    triangle_vertices[3*(area-1) + 1] = j2;
    triangle_vertices[3*(area-1) + 2] = k2;

    area++;
}
void Print(ofstream& file, Grid& grid)
{
    for (int i = 0; i < (3 * grid.triangle); i++)
    {
        file << grid.triangle_vertices[i] << ' ';
        //print 3 vertices per line
        if(i % 3 == 2) file << endl;
    }
}
int Grid::GetMCoordinate(int vertex)
{
    return vertex % m;
}
int Grid::GetNCoordinate(int vertex)
{
    return vertex/m;
}
int Grid::GetN()
{
    return n;
}
int Grid::GetM()
{
    return m;
}
int Grid::GetRectangle()
{
    return rectangle;
}
int Grid::GetTriangle()
{
    return triangle;
}
int Grid::GetNumberofVertices()
{
    return number_of_vertices;
}
