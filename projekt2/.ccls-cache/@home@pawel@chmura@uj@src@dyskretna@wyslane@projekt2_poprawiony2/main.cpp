#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "include/Field.cpp"
#include "include/A2.cpp"


template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os << *ii << " ";
    }
    return os;
}


struct TestResult
{
    bool isHyperbolic;
    int n;
    std::vector<int> line_sizes;
    std::vector<int> possible_counts;
    TestResult() : isHyperbolic(false) {}


    friend std::ostream& operator << (std::ostream& output, const TestResult& res)
    {
        output << res.n << std::endl;
        
        for (int i = 0; i < res.line_sizes.size(); ++i)
        {
            int current_size = res.line_sizes[i];
            for (int j = 0; j < current_size; ++j)
            {
                output << i << " ";
            }
        }
        output << std::endl;
        output << res.possible_counts;
        return output;
    }

    bool operator   ==(const TestResult& res) const {
        return (
            isHyperbolic == res.isHyperbolic &&
            n == res.n &&
            line_sizes == res.line_sizes && 
            possible_counts == res.possible_counts
        );

    }

};


struct LineConfiguration
{
    A2 *plane;
    std::vector<A2::Line*> lines;
    std::vector<A2::Point*> intersections;

    LineConfiguration(A2 *A){ plane = A; }
    bool addLine(A2::Line* line)
    {
        for (A2::Line* l : lines) if (line == l) return false;
        lines.push_back(line);
        
        // compute intersection points
        for (A2::Line* line2 : lines)
        {
            A2::Point* intersection = plane->linesIntersection(line, line2);
            if (intersection != nullptr) addIntersectionPoint(intersection);
        }
        return true;
    }

    void addIntersectionPoint(A2::Point* point)
    {
        for (A2::Point* p : intersections) if (point == p) return;
        intersections.push_back(point);
    }

    // compute signature = (i2, ..., in) - i_k = number of points of intersection crossed by k lines
    std::vector<int> signature() const
    {
        std::vector<int> sig(plane->order(), 0);

        for (A2::Point* intersection : intersections)
        {
            int count = 0;
            for (A2::Line* line : lines)
            {
                if (plane->isPointOnLine(intersection, line)) count++;
            }
            sig[count-2]++;
        }
        return sig;
    }

    friend std::ostream& operator << (std::ostream& output, const LineConfiguration& config)
    {
        output << "Lines: " << std::endl;
        for (A2::Line* line : config.lines)
        {
            output << line->m << ", " << line->c << std::endl; 
        }

        output << "Intersections: " << std::endl;
        for (A2::Point* point : config.intersections)
        {
            output << point << " : " << point->x << ", " << point->y << std::endl;
        }

        std::vector<int> sig = config.signature();
        output << "Signature: " << sig << std::endl;
        output << std::endl;
        return output;
    }
};



int lineSize(A2::Line& line)
{
    int r = 0;
    for (A2::Point* point : line.points)
    {
        if (point->mask == 0) r++;
    }
    return r;
}


int lineIntersectionSize(A2 &plane, A2::Line& line1, A2::Line& line2)
{
    int r = 0;
    std::vector<bool> mem(plane.points.size(), false);
    for (A2::Point* point : line1.points)
    {
        if (point->mask > 0) continue;
        mem[point->id] = true;
    }
    for (A2::Point* point : line2.points)
    {
        if (point->mask > 0) continue;
        if (mem[point->id]) r++;
    }
    return r;
}

/*
Condition 1: by construction
Conditions 2 and 3: this test
Condition 4: implied by construction and conditions 2 and 3 (except if only one line left)
*/
TestResult testConstruction(A2 &plane)
{
    // compute line intersections
    std::vector<std::vector<int>> intersection_matrix(plane.lines.size(), std::vector<int>(plane.lines.size()));
    
    bool failed = false;

    int line_count = 0;

    #pragma omp parallel for
    for (int i = 0; i < plane.lines.size(); ++i)
    {
        A2::Line& line1 = plane.lines[i];
        intersection_matrix[i][i] = lineSize(line1);

        // Condition 3: line not removed and length >= 2 ?
        if (line1.mask == 0)
        {
            if (intersection_matrix[i][i] < 2) failed = true;
            else
            {
                #pragma omp atomic
                line_count++;
            }
        }
        
        for (int j = i+1; j < plane.lines.size(); ++j)
        {
            A2::Line& line2 = plane.lines[j];
            intersection_matrix[i][j] = lineIntersectionSize(plane, line1, line2);
            intersection_matrix[j][i] = intersection_matrix[i][j];
        }
    }
    
    // Condition 4: at lest two lines (of length >= 2)
    if (line_count <= 1 || failed){ return TestResult(); }


    std::vector<bool> possible_line_counts(plane.order()*(plane.order() + 1), false);

    // Condition 2: for each line L, through every point P not on L there are at least two lines disjoint from L
    #pragma omp parallel for
    for (A2::Line& line : plane.lines)
    {
        if (line.mask > 0) continue;    // removed line

        for (A2::Line* parallel : plane.lines_by_slope[line.m])
        {
            if (parallel->mask > 0 || parallel->id == line.id) continue;
            
            for (A2::Point* point_not_on_line : parallel->points)
            {
                if (point_not_on_line->mask > 0) continue;  // removed point

                int lines_through_point_disjoint_from_line_count = 0;

                for (int i = 0; i < plane.order()+1; ++i)
                {
                    A2::Line *line_through_point_not_on_line = point_not_on_line->lines[i];
                    if (line_through_point_not_on_line->mask > 0) continue; // removed line

                    if (intersection_matrix[line.id][line_through_point_not_on_line->id] == 0) lines_through_point_disjoint_from_line_count++;

                }

                // Condition 2
                if (lines_through_point_disjoint_from_line_count < 2)
                {
                    failed = true;
                    //return TestResult();
                }
                possible_line_counts[lines_through_point_disjoint_from_line_count] = true;
            }
        }
    }
    if (failed) return TestResult();

    TestResult r;
    r.isHyperbolic = true;

    // points count
    r.n = 0;
    for (int i = 0; i < plane.points.size(); ++i)
    {
        if (plane.points[i].mask == 0) r.n++;
    }

    // line sizes
    r.line_sizes.resize(plane.order(), 0);
    for (int i = 0; i < intersection_matrix.size(); ++i)
    {
        if (plane.lines[i].mask > 0) continue;
        r.line_sizes[intersection_matrix[i][i]]++;
    }

    // possible line counts
    for (int i = 0; i < possible_line_counts.size(); ++i)
    {
        if (possible_line_counts[i]) r.possible_counts.push_back(i);
    }

 
    return r;
}

std::vector<LineConfiguration> uniqueLineConfigurations(A2 *plane, int k)
{
    if (k == 1)
    {
        LineConfiguration config(plane);
        config.addLine(&plane->lines[0]);
        return std::vector<LineConfiguration>({config});
    }

    std::vector<LineConfiguration> current = uniqueLineConfigurations(plane, k-1);
    std::vector<LineConfiguration> r;
    std::vector<std::vector<int>> signatures;

    for (LineConfiguration& config : current)
    {
        if (config.intersections.size() > 0)
        {
            for (int m = 0; m <= plane->order(); ++m)
            {
                std::vector<bool> cs(plane->order(), false);

                // for each intersection point in this config, note the parameter 'c' of line with clope 'm' through this point
                for (A2::Point* point : config.intersections)
                {
                    cs[point->lines[m]->c] = true;
                }
                // we have eg: cs = 000001000001001001011100011
                // for each interval of zeros in cs, we can add line y=mx+c (c s.t. cs[c] = 0), and the choices are equivalent
                // same for each 1 in cs
                bool add = true;
                std::vector<int> to_add;
                for (int c = 0; c < plane->order(); ++c)
                {
                    if (cs[c] == false)
                    {
                        if (add) to_add.push_back(c);
                        add = false;
                    }
                    if (cs[c] == true)
                    {
                        to_add.push_back(c);
                        add = true;
                    }
                }

                for (int c : to_add)
                {
                    // add line y=mx+c
                    LineConfiguration config_copy = config;
                    if (config_copy.addLine(plane->lines_by_slope[m][c]))
                    {
                        std::vector<int> sig = config_copy.signature();
                        if (std::find(signatures.begin(), signatures.end(), sig) == signatures.end())
                        {
                            //std::cout << "signature " << sig << " not found" << std::endl;
                            r.push_back(config_copy);
                            signatures.push_back(sig);
                        }
                    } 
                }
            }
        }
        else
        {
            // no intersections - all lines parallel
            // case 1 : add a line parallel to these (if possible)
            int m = config.lines[0]->m;

            if (config.lines.size() < plane->order())
            {

                LineConfiguration config_parallel = config;
                std::vector<bool> cs(plane->order(), false);
                for (A2::Line* line : config.lines)
                {
                    cs[line->c] = true;
                }

                for (int i = 0; i < plane->order(); ++i)
                {
                    if (cs[i] == false)
                    {
                        config_parallel.addLine(plane->lines_by_slope[m][i]);
                        break;
                    }
                }
                std::vector<int> config_parallel_sig = config_parallel.signature();
                if (std::find(signatures.begin(), signatures.end(), config_parallel_sig) == signatures.end())
                {
                    r.push_back(config_parallel);
                    signatures.push_back(config_parallel_sig);
                }
            }

            // case 2 : add a non-parallel line
            if (m < plane->order()) m++;
            else m = 0;

            LineConfiguration config_cross = config;
            config_cross.addLine(plane->lines_by_slope[m][0]);
            std::vector<int> config_cross_sig = config_cross.signature();

            if (std::find(signatures.begin(), signatures.end(), config_cross_sig) == signatures.end())
            {
//                std::cout << "signature " << config_cross_sig << " not found" << std::endl;
                r.push_back(config_cross);
                signatures.push_back(config_cross_sig);
            }
        }
    }

    return r;
}


int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cout << "Użycie: " << argv[0] << " typ p a k" << std::endl;
        std::cout << "gdzie: " << std::endl;
        std::cout << "typ = A (przestrzeń afiniczna) lub P (przestrzeń rzutowa)" << std::endl;
        std::cout << "p   - liczba pierwsza" << std::endl;
        std::cout << "a   - liczba naturalna >= 1, wykładnik potęgi p do konstrukcji ciała skończonego F_n gdzie n=p^a" << std::endl;
        std::cout << "k   - liczba prostych do usunięcia przy konstrukcji" << std::endl;
        return 0;
    }

    std::string type = argv[1];
    int p = std::stoi(argv[2]);
    int a = std::stoi(argv[3]);
    int k = std::stoi(argv[4]);
    
    if (k == 0)
    {
        std::cout << "NIE" << std::endl;
    }
    
    if (type == "P")
    {
        if (k == 1)
        {
            // P2 bez jednej prostej to A2 - niehiperboliczna
            std::cout << "NIE" << std::endl;
            return 0;
        }

        // konstrukcja z P dla k prostych = konstrukcja z A dla k-1 prostych 
        k--;
    }
    
    std::vector<TestResult> results;

    A2 A(p,a);
    A.makeLines();

    std::vector<LineConfiguration> configs = uniqueLineConfigurations(&A, k);

    for (LineConfiguration& config : configs)
    {
        for (A2::Line* l : config.lines) A.maskLine(l);
        TestResult r = testConstruction(A);
        for (A2::Line* l : config.lines) A.unmaskLine(l);

        if (!r.isHyperbolic) continue;
        if (results.size() == 0)
        {
            std::cout << "TAK" << std::endl;
        }

        if (std::find_if(results.begin(), results.end(), [r](const TestResult &result)->bool { return result == r; }) == results.end())
        {
            std::cout << r << std::endl;
            results.push_back(r);
        }
    }

    if (results.size() == 0)
    {
        std::cout << "NIE" << std::endl;
    }

    return 0;
}