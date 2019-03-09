#define NONAMELESSUNION
#include <unordered_set>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <filesystem>
#include <iterator>
#include <execution>

#include <ctime>

template <typename Type>
struct Vector : std::vector<Type>
{
    using std::vector<Type>::vector;

    int size() const
    {
        return std::vector<Type>::size();
    }
};

template <typename FirstType, typename SecondType>
using Pair = std::pair<FirstType, SecondType>;

using IntSet = std::unordered_set<int>;
using TagSet = std::unordered_set<std::string>;

template <typename Type>
struct Array2d
{
    int width;
    int height;

    Vector<Type> data;

    struct Row
    {
        Type *begin;
        Type *end;

        Row(Type *begin, Type *end) :
            begin(begin),
            end(end)
        {

        }

        Type &operator [](int idx)
        {
            Type *it = begin + idx;
            if (it < begin | it > end)
            {
                throw std::runtime_error("Column index " + std::to_string(idx) + " out of range");
            }
            return *it;
        }
    };

    struct Iterator //: std::iterator<std::bidirectional_iterator_tag, Row>
    {
        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = int;
        using value_type = int;
        using pointer = value_type * ;
        using reference = value_type & ;

        Array2d<Type> *container;
        int idx;

        Iterator() = default;
        Iterator(Array2d<Type> *container, int idx) :
            container(container),
            idx(idx)
        {

        }

        Row operator *()
        {
            return (*container)[idx];
        }

        Iterator &operator++()
        {
            idx++;
            return *this;
        }

        Iterator &operator--()
        {
            idx--;
            return *this;
        }

        bool operator != (Iterator const &other)
        {
            return idx != other.idx;
        }
    };

    Array2d(int width, int height) :
        width(width),
        height(height)
    {
        data.resize(width * height);
    }

    Row operator [](int idx)
    {
        if (idx < 0 || idx >= height)
        {
            throw std::runtime_error("Row index " + std::to_string(idx) + " out of range");
        }
        return Row(
            data.data() + idx * width,
            data.data() + idx * width + width
        );
    }

    Iterator begin()
    {
        return Iterator(this, 0);
    }

    Iterator end()
    {
        return Iterator(this, height);
    }
};

struct Photo
{
    IntSet tags;
    int idx;
    bool is_vertical;

    Photo() = default;
    Photo(TagSet string_tags, int idx, bool is_vertical) :
        idx(idx),
        is_vertical(is_vertical)
    {
        for (std::string const &tag : string_tags)
        {
            tags.insert(std::hash_value(tag));
        }
    }

    std::string to_string() const
    {
        std::ostringstream oss;
        if (is_vertical)
        {
            oss << 'V';
        }
        else
        {
            oss << 'H';
        }

        oss << ' ' << tags.size();
        for (int tag : tags)
        {
            oss << ' ' << std::to_string(tag);
        }

        return oss.str();
    }
};

using PhotoArray = Vector<Photo>;

struct Slide
{
    PhotoArray photos;
    IntSet tags;

    Slide() = default;
    Slide(PhotoArray photos) :
        photos(std::move(photos))
    {
        for (Photo const &photo : this->photos)
        {
            tags.insert(
                std::begin(photo.tags),
                std::end(photo.tags)
            );
        }
    }
};

using SlideArray = Vector<Slide>;
using IntArray = Vector<int>;

void CombineVerticals(SlideArray &slides, PhotoArray &verticals)
{
    std::sort(
        std::begin(verticals),
        std::end(verticals),
        [](Photo const &lhs, Photo const &rhs)
    {
        return lhs.tags.size() < rhs.tags.size();
    }
    );
    for (int photo_idx = 0; photo_idx < verticals.size() - 1; photo_idx += 2)
    {
        slides.emplace_back(PhotoArray{ verticals[photo_idx], verticals[verticals.size() - photo_idx - 1] });
    }
}

int Score(Slide const &lhs, Slide const &rhs)
{
    int intersection = 0;
    int diff1 = 0;
    int diff2 = 0;

    for (int tag : lhs.tags)
    {
        if (rhs.tags.count(tag))
        {
            intersection++;
        }
        else
        {
            diff1++;
        }
    }

    for (int tag : rhs.tags)
    {
        if (!lhs.tags.count(tag))
        {
            diff2++;
        }
    }

    return std::min<int>(
        {
            intersection,
            diff1,
            diff2
        }
    );
}

int Score(SlideArray const &slides)
{
    int sum = 0;
    for (int idx = 0; idx < slides.size() - 1; idx++)
    {
        sum += Score(
            slides[idx],
            slides[idx + 1]
        );
    }
    return sum;
}

SlideArray SortSlides(SlideArray const &slides)
{
    IntArray output_indices;
    output_indices.reserve(slides.size());

    std::srand(std::time(0));

    int first_idx = std::rand() % slides.size();
    output_indices.push_back(first_idx);

    IntSet indices;
    for (int idx = 0; idx < slides.size(); idx++)
    {
        if (first_idx != idx)
        {
            indices.insert(idx);
        }
    }

    while (output_indices.size() < slides.size())
    {
        int best_score = -1;
        int best_idx = -1;

        for (int idx : indices)
        {
            int s = Score(
                slides[output_indices.back()],
                slides[idx]
            );

            if (best_idx == -1 || s > best_score)
            {
                best_score = s;
                best_idx = idx;
            }
        }

        if (indices.count(best_idx))
        {
            indices.erase(best_idx);
            output_indices.push_back(best_idx);
        }

        if (output_indices.size() > 2)
        {
            while (Score(
                slides[output_indices[output_indices.size() - 2]],
                slides[output_indices[output_indices.size() - 1]]
            ) < Score(
                slides[output_indices[output_indices.size() - 1]],
                slides[output_indices[0]]
            ))
            {
                output_indices.insert(
                    std::begin(output_indices),
                    output_indices[output_indices.size() - 1]
                );
                output_indices.pop_back();
            }
        }
    }

    SlideArray output;
    output.reserve(slides.size());

    for (int idx : output_indices)
    {
        output.push_back(slides[idx]);
    }

    return output;
}

void SolveProblem(std::filesystem::path filename, int nthreads = 1)
{
    PhotoArray photos;

    std::fstream file(filename);

    int photo_count;
    file >> photo_count;

    for (int slide_idx = 0; slide_idx < photo_count; slide_idx++)
    {
        std::string v;
        file >> v;

        int tag_count;
        file >> tag_count;

        TagSet set;
        while (tag_count--)
        {
            std::string tag;
            file >> tag;
            set.insert(tag);
        }

        photos.emplace_back(set, slide_idx, v == "V");
    }

    PhotoArray verticals;
    SlideArray slides;

    for (Photo const &photo : photos)
    {
        if (photo.is_vertical)
        {
            verticals.push_back(photo);
        }
        else
        {
            slides.emplace_back(PhotoArray{ photo });
        }
    }

    CombineVerticals(slides, verticals);

    std::sort(
        std::begin(slides),
        std::end(slides),
        [](Slide const &lhs, Slide const &rhs)
    {
        return lhs.tags.size() < rhs.tags.size();
    }
    );

    //Array2d<int> order(slides.size() - 1, slides.size());
    //Array2d<int> scores(slides.size(), slides.size());
    //std::for_each(
    //    std::execution::par_unseq,
    //    std::begin(order),
    //    std::end(order),
    //    [&](Array2d<int>::Row row)
    //{
    //    int row_idx = std::distance(order.data.data(), row.begin) / order.width;

    //    int *it = row.begin;

    //    Slide const &slide = slides[row_idx];

    //    int slide_idx = slides.size();
    //    while (slide_idx--)
    //    {
    //        if (slide_idx != row_idx)
    //        {
    //            *it++ = slide_idx;
    //            scores[row_idx][slide_idx] = Score(slide, slides[slide_idx]);
    //        }
    //        else
    //        {
    //            scores[row_idx][slide_idx] = 0;
    //        }
    //    }

    //    std::sort(row.begin, row.end, [row_idx, &scores](int lhs, int rhs)
    //    {
    //        return scores[row_idx][lhs] < scores[row_idx][lhs];
    //    });
    //});

    Vector<Pair<int, SlideArray>> splits;

    int step = slides.size() / nthreads;
    for (int idx = 0; idx < nthreads; idx++)
    {
        int start = step * idx;
        int end = step * (idx + 1);
        if (idx == nthreads - 1)
        {
            end = slides.size();
        }
        splits.emplace_back(
            idx,
            SlideArray(
                std::next(std::begin(slides), start),
                std::next(std::begin(slides), end)
            )
        );
    }

    std::for_each(
        std::execution::par_unseq,
        std::begin(splits),
        std::end(splits),
        [step, &slides](Pair<int, SlideArray> const &pair)
    {
        SlideArray sorted_slides = SortSlides(pair.second);
        std::copy(
            std::begin(sorted_slides),
            std::end(sorted_slides),
            std::next(
                std::begin(slides),
                step * pair.first
            )
        );
    }
    );

    int score = Score(slides);
    std::cout << score << std::endl;

    std::filesystem::path dir = filename.stem();
    create_directories(dir);

    std::ofstream outfile(dir / (std::to_string(score) + ".txt"));

    outfile << slides.size() << std::endl;
    for (Slide const &slide : slides)
    {
        for (Photo const &photo : slide.photos)
        {
            outfile << photo.idx << ' ';
        }
        outfile << std::endl;
    }
}

int main(int argc, char *argv[])
{
    //SolveProblem("b_lovely_landscapes.txt", 1);
    //SolveProblem("e_shiny_selfies.txt", 1);
    //SolveProblem("d_pet_pictures.txt", 1);
    //SolveProblem("c_memorable_moments.txt", 1);

    for (int idx = 1; idx < argc; idx++)
    {
        switch (argv[idx][0])
        {
        case 'a':
        case 'A':
            SolveProblem("a_example.txt");
            break;
        case 'b':
        case 'B':
            SolveProblem("b_lovely_landscapes.txt");
            break;
        case 'c':
        case 'C':
            SolveProblem("c_memorable_moments.txt");
            break;
        case 'd':
        case 'D':
            SolveProblem("d_pet_pictures.txt");
            break;
        case 'e':
        case 'E':
            SolveProblem("e_shiny_selfies.txt");
        default:
            break;
        }
    }

    return 0;
}