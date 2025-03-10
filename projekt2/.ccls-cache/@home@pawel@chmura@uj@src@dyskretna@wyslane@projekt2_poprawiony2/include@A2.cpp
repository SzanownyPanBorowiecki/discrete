#ifndef A2_CPP
#define A2_CPP

#include <vector>
#include <utility>
#include "Field.cpp"

struct A2
{
    struct Line;

    struct Point
    {
        Field x; Field y;
        int id;
        int mask;

        std::vector<Line*> lines;
        Point() : 
            lines(Field::size() + 1),
            mask(0)
            {};

        Point(Field cx, Field cy) : 
            lines(Field::size() + 1),
            mask(0)
            { x=cx; y=cy; }

        operator int() const { return Field::size()*int(y)+int(x); }
        bool operator   ==(const Point& p) const { return x == p.x && y == p.y; }

        
        friend std::ostream& operator << (std::ostream& os, const Point& p)
        {
            os << "(" << p.x << ", " <<  p.y << ")";
            return os;
        }
    };

    struct Line
    {
        int id;
        int m, c;
        int mask;

        std::vector<Point*> points;
        Line() : points(Field::size()), mask(0) {}
    };

    int ord;
    std::vector<Point> points;
    std::vector<Line> lines;

    std::vector<std::vector<Line*>> lines_by_slope;
    
    // multiplication and addition tables for speed
    std::vector<std::vector<int>> mult;
    std::vector<std::vector<int>> add;


    void makeLines()
    {
        int n = Field::size();
        
        // non-vertical
        #pragma omp parallel for
        for (int m = 0; m < n; ++m)
        {
            for (int c = 0; c < n; ++c)
            {
                int lineid = m*n + c;
                lines[lineid].id = lineid;
                lines[lineid].m = m;
                lines[lineid].c = c;
                
                for (int x = 0; x < n; ++x)
                {
                    int y = add[mult[m][x]][c];
                    int pointid = n*y + x;
                    lines[lineid].points[x] = &points[pointid];
                    points[pointid].lines[m] = &lines[lineid];
                }
                lines_by_slope[m][c] = &lines[lineid];

            }
        }

        // vertical
        #pragma omp parallel for
        for (int x = 0; x < n; ++x)
        {
            int lineid = n*n + x;
            lines[lineid].id = lineid;
            lines[lineid].m = n;
            lines[lineid].c = x;

            // add points (x, {0,1,2,....})
            for (int y = 0; y < n; ++y)
            {
                int pointid = n*y + x;
                lines[lineid].points[y] = &points[pointid];
                points[pointid].lines[n] = &lines[lineid];
            }

            lines_by_slope[n][x] = &lines[lineid];
        }

    }

    A2(int p, int e)
    {
        Field::init(p,e);
        int n = Field::size();
        ord = n;
        points.resize(n*n);
        lines.resize(n*n + n);
        lines_by_slope.resize(n+1, std::vector<Line*>(n));
        
        mult.resize(n, std::vector<int>(n, 0));
        add.resize(n, std::vector<int>(n, 0));

        //std::cout << "make multiplication tables" << std::endl;
        for (int i = 0; i < n; ++i)
        {
            Field x(i);
            #pragma omp parallel for
            for (int j = 0; j < n; ++j)
            {
                Field y(j);
                if (j >= i)
                {
                    mult[i][j] = int(x*y);
                    add[i][j] = int(x+y);
                }
                else
                {
                    mult[i][j] = mult[j][i];
                    add[i][j] = add[j][i];
                }

                int pointid = j*Field::size() + i; 
                points[pointid] = Point(x, y);
                points[pointid].id = pointid;
            }
        }
        //std::cout << "done A2 init" << std::endl;
    }

    int order(){ return Field::size(); }

    void maskLine(Line* line)
    {
        if (line->mask > 0) return;
        line->mask++;
        // increase mask level for each point (more mask -> removed by more lines)
        for (A2::Point* point : line->points)
        {
            point->mask++;
        }
    }

    void unmaskLine(Line* line)
    {
        if (line->mask == 0) return;
        line->mask--;
        for (A2::Point* point : line->points)
        {
            if (point->mask > 0) point->mask--;
        }
    }

    Point* linesIntersection(Line* line1, Line* line2)
    {
        //std::cout << "linesIntersection(" << line1 << ", " << line2 << ")" << std::endl;
        if (line1->m == line2->m) return nullptr;
        if (line1->m == order())
        {
            int x = line1->c;
            int y = add[mult[line2->m][x]][line2->c];
            int pointid = order()*y + x;
            //intersection = ;
            //std::cout << "intersection = " << intersection << std::endl;
            return &points[pointid];
        }
        else if (line2->m == order()) return linesIntersection(line2, line1);
        else
        {
            Field c1(line1->c);
            Field c2(line2->c);
            Field m1(line1->m);
            Field m2(line2->m);
            Field fx = (c2-c1)/(m1-m2);
            Field fy = m1*fx + c1;
            int x = int(fx);
            int y = int(fy);
            int pointid = order()*y + x;

            //intersection = &points[pointid];
            //std::cout << "intersection = " << intersection << std::endl;

            return &points[pointid];
        }
    }

    bool isPointOnLine(Point* point, Line* line)
    {
        int px = int(point->x);
        int py = int(point->y);
        if (line->m == order()) return (px == line->c);
        return py == add[mult[line->m][px]][line->c];
    }
};

#endif