#include "utils.h"
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#define open _open
#define close _close
#define O_WRONLY _O_WRONLY
#define O_CREAT _O_CREAT
#define O_EXCL _O_EXCL
#define S_IRUSR _S_IREAD
#define S_IWUSR _S_IWRITE
#else
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#endif

int openWithFcntl(const char* filename) {
    // Use the global scope "::open()" to specify fcntl version
    int fd = ::open(filename, O_WRONLY | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);

    if (fd == -1) {
        std::cerr << "Error: Could not open file with fcntl\n";
    }
    else {
        std::cout << "File opened successfully with fcntl, fd: " << fd << "\n";
        close(fd); // Always close file descriptors when done
    }
    return fd;
}

string replace_char(string s, char c1, char c2) {
    int l = s.length();
    // loop to traverse in the string
    for (int i = 0; i < l; i++) {
        // check for c1 and replace
        if (s[i] == c1)
            s[i] = c2;

        // check for c2 and replace
        else if (s[i] == c2)
            s[i] = c1;
    }
    return s;
}

vector<string> split(string str, const string delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
    str.erase(remove(str.begin(), str.end(), '\r'), str.end());
    str.erase(remove(str.begin(), str.end(), '\n'), str.end());

    do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}

vector<string> split(string str)
{
	istringstream iss(str);
	vector<string> results((istream_iterator<string>(iss)), istream_iterator<string>());
	return results;
}

vector<int> vect_add(vector<int>& vect1, vector<int>& vect2) {
    if (vect1.size() != vect2.size()) {
		cout << "Error: vectors must have the same size" << endl;
		exit(1);
	}
    vector<int> temp;
    for (int i = 0; i < vect1.size(); i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

vector<float> vect_subtract(vector<float>& vect1, vector<float>& vect2) {
    if (vect1.size() != vect2.size()) {
        cout << "Error: vectors must have the same size" << endl;
        exit(1);
    }
    vector<float> temp;
    for (int i = 0; i < vect1.size(); i++) {
        temp.push_back(vect1[i] - vect2[i]);
    }
    return temp;
}

vector<float> vect_add(vector<float>& vect1, vector<float>& vect2){
    if (vect1.size() != vect2.size()) {
        cout << "Error: vectors must have the same size" << endl;
        exit(1);
    }
    vector<float> temp;
    for (int i = 0; i < vect1.size(); i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

vector<double> vect_add(vector<double>& vect1, vector<double>& vect2){
    if (vect1.size() != vect2.size()) {
        cout << "Error: vectors must have the same size" << endl;
        exit(1);
    }
    vector<double> temp;
    for (int i = 0; i < vect1.size(); i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

vector<double> vect_add(vector<double>& vect1, vector<float>& vect2){
    if (vect1.size() != vect2.size()) {
        cout << "Error: vectors must have the same size" << endl;
        exit(1);
    }
    vector<double> temp;
    for (int i = 0; i < vect1.size(); i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

vector<float> vect_add(vector<float>& vect1, float vect2[3]) {
    if (vect1.size() != 3) {
        cout << "Error: vectors must have the same size" << endl;
        exit(1);
    }
    vector<float> temp;
    for (int i = 0; i < 3; i++) {
        temp.push_back(vect1[i] + vect2[i]);
    }
    return temp;
}

int get_index(vector<int>& vect, int elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<string>& vect, string elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int get_index(vector<float>& vect, float elem) {
    auto it = find(vect.begin(), vect.end(), elem);
    if (it != vect.end()) { // If element was found calc the index of K
        int index = it - vect.begin();
        return index;
    }
    else { // If the element is not present in the vector
        return -1;
    }
}

int vect_max(vector<int>& vect) {
    int max = vect[0];
    for (int i : vect) { if (i > max) { max = i; } }
    return max;
}

float vect_mat(vector<float>& vect) {
    float max = vect[0];
    for (float i : vect) { if (i > max) { max = i; } }
    return max;
}

int kron_del(int int1, int int2) {
    if (int1 == int2) { return 1; }
    else { return 0; }
}

int kron_del(float float1, float float2) {
    if (round_to(float1, 5) == round_to(float2, 5)) { return 1; }
    else { return 0; }
}

int kron_del(int int1, float float2){
    if (float(int1) == round_to(float2, 5)) { return 1; }
    else { return 0; }
}

int kron_del(float float1, int int2){
    if (round_to(float1, 5) == float(int2)) { return 1; }
    else { return 0; }
}

float round_to(float val, int digits) {
    float mult = pow(10.0, digits);
    return roundf(val * mult) / mult;
}

float vect_max(vector<float>& vect) {
    float max = vect[0];
    for (float i : vect) { if (i > max) { max = i; } }
    return max;
}

vector<int> vect_permut(vector<int>& vect) {
    vector<int> vect_temp;
    vect_temp.assign(vect.begin(), vect.end());
    vector<pair<int, int>> indices;
    for (int i = 0; i < vect_temp.size(); ++i) {
        indices.push_back(make_pair(vect_temp[i], i));
    }
    sort(indices.begin(), indices.end());

    vector<int> permut(vect.size());
    for (int i = 0; i < indices.size(); ++i) {
        permut[i] = indices[i].second;
    }
    return permut;
}

vector<int> vect_permut(vector<float>& vect) {
    vector<int> vect_temp;
    vect_temp.assign(vect.begin(), vect.end());
    vector<pair<int, int>> indices;
    for (int i = 0; i < vect_temp.size(); ++i) {
        indices.push_back(make_pair(vect_temp[i], i));
    }
    sort(indices.begin(), indices.end());

    vector<int> permut(vect.size());
    for (int i = 0; i < indices.size(); ++i) {
        permut[i] = indices[i].second;
    }
    return permut;
}

void sort_vect(vector<int>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<int> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<float>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<float> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<vector<float>>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<vector<float>> temp;
    temp = vect;// .insert(temp.end(), vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

void sort_vect(vector<vector<int>>& vect, vector<int>& perm) {
    // Rearrange the vector based on the sorting permutation
    vector<vector<int>> temp;
    temp.assign(vect.begin(), vect.end());
    for (int i = 0; i < vect.size(); ++i) {
        vect[i] = temp[perm[i]];
    }
}

vector<float> pos_transform(vector<float>& pos, vector<vector<float>>& trans) {
    //cartesian coordinate vector x*a+y*b+z*c
    vector<float> vect(pos.size());
    vect[0] += pos[0] * trans[0][0] + pos[1] * trans[1][0] + pos[2] * trans[2][0];
    vect[1] += pos[0] * trans[0][1] + pos[1] * trans[1][1] + pos[2] * trans[2][1];
    vect[2] += pos[0] * trans[0][2] + pos[1] * trans[1][2] + pos[2] * trans[2][2];
    return vect;
}

vector<float> scale_vect(vector<float>& vect, float numb) {
    vector<float> new_vect(vect.size(), 0.0);
    for (int i = 0; i < vect.size(); i++) { new_vect[i] = numb * vect[i]; }
    return new_vect;
}

vector<int> scale_vect(vector<int>& vect, int numb) {
    vector<int> new_vect(vect.size(), 0);
    for (int i = 0; i < vect.size(); i++) { new_vect[i] = numb * vect[i]; }
    return new_vect;
}

int sgn(float v) {
    return (v > 0) - (v < 0);
}

int sign(float v) {
    if (v >= 0) { return 1; }
    else { return -1; }
}

float dot_pos(vector<float>& vect1, vector<float>& vect2) {
    float prod = 0;
    for (int i = 0; i < 3; i++) {
        prod += vect1[i] * vect2[1];
    }
    return prod;
}

float pos_dist(vector<float> pos1, vector<float> pos2) {
    float dist = 0.0;
    for (int i = 0; i < 3; i++) { dist += pow(pos1[i] - pos2[i], 2); }
    return sqrt(dist);
}

vector<float> pos_round(vector<float>& pos, int digits) {
	vector<float> new_pos;
	for (int i = 0; i < 3; i++) { new_pos.push_back(round_to(pos[i], digits)); }
	return new_pos;
}

bool fcomp(float a, float b, float tol) {
	return (fabs(a - b) < tol);
}

bool pos_comp(vector<float>& pos1, vector<float>& pos2, float tol) {
    for (int i = 0; i < 3; i++) {
		if (!fcomp(pos1[i], pos2[i], tol)) { return false; }
	}
	return true;
}

vector<float> cart_to_frac(vector<float>& pos, vector<vector<float>>& lat_vec) {
    vector<float> new_pos {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++){
            new_pos[i] += pos[i] * lat_vec[j][i];
        }
    }
    return new_pos;
}
 
vector<float> frac_to_cart(vector<float>& pos, vector<vector<float>>& lat_vec) {
    vector<float> new_pos {0, 0, 0};
    // unfinished
    return new_pos;
}

//vector<float> matrix_inv(Matrix33d m) {
//    // computes the inverse of a matrix m
//    double det = m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
//                 m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
//                 m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
//
//    double invdet = 1 / det;
//
//    Matrix33d minv; // inverse of matrix m
//    minv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * invdet;
//    minv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invdet;
//    minv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invdet;
//    minv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invdet;
//    minv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invdet;
//    minv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invdet;
//    minv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * invdet;
//    minv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invdet;
//    minv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invdet;
//    
//}
