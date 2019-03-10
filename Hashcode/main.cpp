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
#include <random>

#include <ctime>

template <typename Type>
struct Vector : std::vector<Type>
{
	using std::vector<Type>::vector;

	long long size() const
	{
		return (long long)std::vector<Type>::size();
	}
};

template <typename FirstType, typename SecondType>
using Pair = std::pair<FirstType, SecondType>;

using IntSet = std::unordered_set<int>;
using LongLongSet = std::unordered_set<long long>;
using TagSet = std::unordered_set<std::string>;

template <typename Type>
struct Array2d
{
	int width;
	int height;

	Vector<Type> data;

	struct Row
	{
		Type *from;
		Type *to;

		Row() = default;
		Row(Type *begin, Type *end) :
			from(begin),
			to(end)
		{
			if (end < begin)
			{
				throw std::runtime_error("End if before begin");
			}
		}

		Type &operator [](int idx)
		{
			Type *it = from + idx;
			if (it < from || it >= to)
			{
				throw std::runtime_error("Column index " + std::to_string(idx) + " out of range");
			}
			return *it;
		}

		Type const &operator [](int idx) const
		{
			Type *it = from + idx;
			if (it < from || it >= to)
			{
				throw std::runtime_error("Column index " + std::to_string(idx) + " out of range");
			}
			return *it;
		}

		Type *begin()
		{
			return from;
		}

		Type const *begin() const
		{
			return from;
		}

		Type *end()
		{
			return to;
		}

		Type const *end() const
		{
			return to;
		}
	};

	struct Iterator //: std::iterator<std::bidirectional_iterator_tag, Row>
	{
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type = long long;
		using value_type = int;
		using pointer = value_type * ;
		using reference = value_type & ;

		Array2d<Type> *container;
		int idx;

		Iterator() = default;
		Iterator(Array2d<Type> *container, long long idx) :
			container(container),
			idx(idx)
		{

		}

		Row operator *() const
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

	Array2d(int width, int height)
	{
		resize(width, height);
	}

	void resize(int new_width, int new_height)
	{
		width = new_width;
		height = new_height;
		data.resize(new_width * (long long)new_height, Type{});
	}

	Row operator [](int idx)
	{
		if (idx < 0 || idx >= height)
		{
			throw std::runtime_error("Row index " + std::to_string(idx) + " out of range");
		}
		return Row(
			data.data() + idx * (long long)width,
			data.data() + (idx + 1) * (long long)width
		);
	}

	void clear()
	{
		data.clear();
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
	LongLongSet tags;
	int idx;
	bool is_vertical;

	Photo() = default;
	Photo(TagSet string_tags, int idx, bool is_vertical) :
		idx(idx),
		is_vertical(is_vertical)
	{
		for (std::string const &tag : string_tags)
		{
			tags.insert((long long)std::hash_value(tag));
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
		for (long long tag : tags)
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
	LongLongSet tags;

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
using DoubleArray = Vector<double>;

int CommonTagCount(Photo const &lhs, Photo const &rhs)
{
	int count = 0;
	for (int tag : lhs.tags)
	{
		if (rhs.tags.count(tag))
		{
			count++;
		}
	}
	return count;
}

Array2d<int> order(1000, 1000);
Array2d<int> scores(1000, 1000);

void CombineVerticals(SlideArray &slides, PhotoArray &photos)
{
	std::shuffle(
		std::begin(photos),
		std::end(photos),
		std::default_random_engine()
	);

	order.resize(photos.size() - 1, photos.size());
	scores.resize(photos.size(), photos.size());

	IntArray score_sums;
	score_sums.resize(photos.size());

	DoubleArray score_averages;
	score_averages.resize(photos.size());

	IntSet used_photos;

	std::cout << "Precalculating vertical photo scores" << std::endl;
	std::for_each(
		std::execution::par_unseq,
		std::begin(order),
		std::end(order),
		[&](Array2d<int>::Row row)
	{
		int row_idx = (int)(std::abs(order.data.data() - row.from) / order.width);

		Photo const &photo = photos[row_idx];

		int column_idx = (int)photos.size() - 1;
		while (column_idx--)
		{
			order[row_idx][column_idx] = column_idx;
			if (column_idx >= row_idx)
			{
				order[row_idx][column_idx] += 1;
			}
		}

		int photo_idx = (int)photos.size();
		while (photo_idx--)
		{
			int score = CommonTagCount(photo, photos[photo_idx]);

			scores[row_idx][photo_idx] = score;
		}
		scores[row_idx][row_idx] = 0;

		score_sums[row_idx] = std::accumulate(
			std::begin(scores[row_idx]),
			std::end(scores[row_idx]),
			0
		);

		score_averages[row_idx] = score_sums[row_idx] / (double)order.width;

		std::sort(row.from, row.to, [row_idx, &score_averages](int lhs, int rhs)
		{
			double dlhs = std::abs(scores[row_idx][lhs]);
			double drhs = std::abs(scores[row_idx][rhs]);
			return dlhs > drhs;
		});

		//for (int idx = 0; idx < order.width; idx += 2)
		//{
		//	std::swap(row[idx], row[order.width - idx - 1]);
		//}
	});

	std::cout << "Combining " << photos.size() << " vertical photos" << std::endl;
	for (int first_idx = 0; first_idx < order.height; first_idx++)
	{
		if (used_photos.count(first_idx))
		{
			continue;
		}

		int second_idx = -1;

		for (int column_idx = 0; column_idx < order.width; column_idx++)
		{
			int photo_idx = order[first_idx][column_idx];
			if (first_idx == photo_idx)
			{
				continue;
			}

			if (used_photos.count(photo_idx))
			{
				continue;
			}

			second_idx = photo_idx;
		}

		if (second_idx != -1)
		{
			slides.emplace_back(PhotoArray{ photos[first_idx], photos[second_idx] });

			used_photos.insert(first_idx);
			used_photos.insert(second_idx);
		}
	}
}

int Score(Slide const &lhs, Slide const &rhs)
{
	int intersection = 0;
	int diff1 = 0;
	int diff2 = 0;

	for (long long tag : lhs.tags)
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

	for (long long tag : rhs.tags)
	{
		if (!lhs.tags.count(tag))
		{
			diff2++;
		}
	}

	return std::min(
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

int Score(SlideArray const &slides, Array2d<int>::Row row)
{
	int sum = 0;
	for (int idx = 0; idx < slides.size() - 1; idx++)
	{
		sum += Score(
			slides[row[idx]],
			slides[row[idx + 1]]
		);
	}
	return sum;
}

void FindBest(
	Array2d<int> &order,
	Array2d<int> &scores,
	IntSet const &potential_indices,
	int first_idx,
	int skip,
	int depth,
	int max_depth,
	int total_score,
	IntArray &indices,
	int &best_score,
	Array2d<int>::Row best_combination
)
{
	int count = 0;
	bool wraparound = false;
	for (int column_idx = skip; count != 3 && count != potential_indices.size(); column_idx++)
	{
		if (column_idx >= order.width)
		{
			column_idx -= order.width;
			wraparound = true;
		}

		if (wraparound && skip == column_idx)
		{
			break;
		}

		int slide_idx = order[first_idx][column_idx];

		auto begin = std::begin(indices);
		auto end = std::next(std::begin(indices), depth);

		int score = scores[first_idx][slide_idx];

		bool not_used = std::find(begin, end, slide_idx) == end;
		if (potential_indices.count(slide_idx) && (not_used || potential_indices.size() == 1))
		{
			indices[depth] = slide_idx;
			count++;

			if (depth < max_depth)
			{
				FindBest(
					order,
					scores,
					potential_indices,
					indices[depth],
					skip,
					depth + 1,
					max_depth,
					total_score + score,
					indices,
					best_score,
					best_combination
				);
			}
			else if (score + total_score > best_score)
			{
				best_score = score + total_score;

				std::copy(
					begin,
					end,
					std::begin(best_combination)
				);
			}
		}
	}
}

void FindBestCombination(
	Array2d<int> &order,
	Array2d<int> &scores,
	IntSet const &potential_indices,
	int first_idx,
	int depth,
	int max_depth,
	int &best_score,
	Vector<int> &best_indices
)
{
	Vector<int> candidates;
	candidates.reserve(max_depth);

	for (int column_idx = 0; column_idx < order.width && candidates.size() < depth; column_idx++)
	{
		int slide_idx = order[first_idx][column_idx];
		if (potential_indices.count(slide_idx))
		{
			candidates.push_back(slide_idx);
		}
	}

	for (int start_idx : candidates)
	{
		Vector<int> next_n;
		next_n.reserve(max_depth);

		next_n.push_back(start_idx);

		int score = scores[first_idx][start_idx];

		int last_idx = start_idx;
		while (next_n.size() < depth - 1)
		{
			for (int column_idx = 0; column_idx < order.width; column_idx++)
			{
				int slide_idx = order[last_idx][column_idx];
				if (potential_indices.count(slide_idx) && std::find(std::begin(next_n), std::end(next_n), slide_idx) == std::end(next_n))
				{
					next_n.push_back(slide_idx);
					score += scores[last_idx][slide_idx];
					last_idx = slide_idx;
					break;
				}
			}
		}

		if (best_score == -1 || best_score < score)
		{
			best_score = score;
			best_indices = next_n;
		}
	}
}

std::random_device rd;std::mt19937 engine(rd());

void SortSlides(SlideArray const &slides, Array2d<int> &scores, Array2d<int> &order, Array2d<int>::Row dest)
{
	int first_idx = std::uniform_int_distribution<int>(0, slides.size())(engine);
	int last_idx = first_idx;

	IntSet potential_indices;
	for (int idx = 0; idx < slides.size(); idx++)
	{
		if (last_idx != idx)
		{
			potential_indices.insert(idx);
		}
	}

	//int skip = 0;
	int max_depth = 1;
	//int last_size = 0;

	Array2d<int> best_matrix(max_depth + 1, potential_indices.size());
	IntArray best_vector;

	int *it = dest.from;
	*it++ = last_idx;
	while (it != dest.to)
	{
		int best_back_idx = -1;
		int best_back_score = -1;
		for (int column_idx = 0; column_idx < order.width; column_idx++)
		{
			int slide_idx = order[last_idx][column_idx];
			if (potential_indices.count(slide_idx))
			{
				best_back_idx = slide_idx;
				best_back_score = scores[last_idx][slide_idx];
				break;
			}
		}

		int best_front_idx = -1;
		int best_front_score = -1;
		for (int column_idx = 0; column_idx < order.width; column_idx++)
		{
			int slide_idx = order[first_idx][column_idx];
			if (potential_indices.count(slide_idx))
			{
				best_front_idx = slide_idx;
				best_front_score = scores[first_idx][slide_idx];
				break;
			}
		}

		if (best_back_score > best_front_score)
		{
			if (potential_indices.count(best_back_idx))
			{
				potential_indices.erase(best_back_idx);
				*it++ = best_back_idx;
				last_idx = best_back_idx;
			}
		}
		else
		{
			if (potential_indices.count(best_front_idx))
			{
				potential_indices.erase(best_front_idx);
				std::copy(dest.from, it, dest.from + 1);
				it++;
				dest[0] = best_front_idx;
				first_idx = best_front_idx;
			}
		}

		//if (std::abs(it - dest.from) > 2)
		//{
		//	int second_last_idx = *(it - 2);
		//	while (scores[second_last_idx][last_idx] < scores[last_idx][first_idx])
		//	{
		//		std::copy(dest.from, it - 1, dest.from + 1);
		//		dest[0] = last_idx;
		//		first_idx = last_idx;
		//		last_idx = second_last_idx;
		//	}
		//}
	}
}

void SolveProblem(std::filesystem::path filename, int ntimes = 1)
{
	std::fstream file(filename);

	file.seekg(0, file.end);
	size_t file_size = file.tellg();
	file.seekg(0, file.beg);

	std::string input;
	input.resize(file_size);

	file.read(&input[0], file_size);

	std::istringstream iss(input);

	int photo_count;
	iss >> photo_count;

	PhotoArray verticals;
	verticals.reserve(photo_count / 2);

	SlideArray slides;
	slides.reserve(photo_count / 2);

	for (int photo_idx = 0; photo_idx < photo_count; photo_idx++)
	{
		std::string v;
		iss >> v;

		int tag_count;
		iss >> tag_count;

		TagSet set;
		while (tag_count--)
		{
			std::string tag;
			iss >> tag;
			set.insert(tag);
		}

		Photo photo(set, photo_idx, v == "V");

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

	order.resize(slides.size() - 1, slides.size());
	scores.resize(slides.size(), slides.size());

	std::cout << "Precalculating slide scores" << std::endl;

	std::for_each(
		std::execution::par_unseq,
		std::begin(order),
		std::end(order),
		[&](Array2d<int>::Row row)
	{
		int row_idx = (int)(std::abs(order.data.data() - row.from) / order.width);

		Slide const &slide = slides[row_idx];

		int column_idx = (int)slides.size() - 1;
		while (column_idx--)
		{
			order[row_idx][column_idx] = column_idx;
			if (column_idx >= row_idx)
			{
				order[row_idx][column_idx] += 1;
			}
		}

		int slide_idx = (int)slides.size();
		while (slide_idx--)
		{
			int score = Score(slide, slides[slide_idx]);

			scores[row_idx][slide_idx] = score;
		}
		scores[row_idx][row_idx] = 0;

		int score_sum = std::accumulate(
			std::begin(scores[row_idx]),
			std::end(scores[row_idx]),
			0
		);

		int score_max = *std::max_element(
			std::begin(scores[row_idx]),
			std::end(scores[row_idx])
		);

		double score_average = score_sum / (double)order.width;

		std::sort(
			std::begin(row),
			std::end(row),
			[row_idx, score_average](int lhs, int rhs)
		{
			double dlhs = std::abs(scores[row_idx][lhs]);
			double drhs = std::abs(scores[row_idx][rhs]);
			return dlhs > drhs;
		});

		//std::shuffle(
		//	std::begin(row),
		//	std::end(row),
		//	std::default_random_engine()
		//);

		//for (int idx = 0; idx < order.width; idx += 2)
		//{
		//	std::swap(row[idx], row[order.width - idx - 1]);
		//}
	});

	std::cout << "Sorting " << slides.size() << " slides" << std::endl;

	Array2d<int> outputs(slides.size(), ntimes);
	IntArray output_scores(ntimes);
	std::for_each(
		std::execution::par_unseq,
		std::begin(outputs),
		std::end(outputs),
		[&](Array2d<int>::Row row)
	{
		int row_idx = (int)(std::abs(outputs.data.data() - row.from) / outputs.width);
		SortSlides(slides, scores, order, row);
		output_scores[row_idx] = Score(slides, row);
	});

	for (int idx = 0; idx < ntimes; idx++)
	{
		Array2d<int>::Row row = outputs[idx];
		int score = output_scores[idx];

		std::cout << score << std::endl;

		std::filesystem::path dir = filename.stem();
		create_directories(dir);

		std::ofstream outfile(dir / (std::to_string(score) + ".txt"));

		outfile << slides.size() << std::endl;
		for (int slide_idx : row)
		{
			for (Photo const &photo : slides[slide_idx].photos)
			{
				outfile << photo.idx << ' ';
			}
			outfile << std::endl;
		}
	}
}

int main(int argc, char *argv[])
try
{
	std::srand((unsigned)std::time(0));
	for (int idx = 1; idx < argc; idx++)
	{
		switch (argv[idx][0])
		{
		case 'a':
		case 'A':
			SolveProblem("a_example.txt", 100);
			break;
		case 'b':
		case 'B':
			SolveProblem("b_lovely_landscapes.txt", 100);
			break;
		case 'c':
		case 'C':
			SolveProblem("c_memorable_moments.txt", 100);
			break;
		case 'd':
		case 'D':
			SolveProblem("d_pet_pictures.txt", 100);
			break;
		case 'e':
		case 'E':
			SolveProblem("e_shiny_selfies.txt", 100);
		default:
			break;
		}
	}

	return 0;
}
catch (std::exception const &e)
{
	std::cout << e.what() << std::endl;
	return 1;
}