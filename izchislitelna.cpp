#include <iostream>
#include <sstream>
#include <string>
#include <hash_map>
#include <cstring>
#include <vector>
#include <algorithm>

using namespace std;

#define MAX_SIZE 70
#define LAT_SIZE 140

/* This matrix is used in the sequence alignment algorythm.
 * if the residue at position i of sequence #1 is the same as the residue at position j of sequence #2: sim_mat[i][j] = 2 (match score)
 * otherwise: sim_mat[i][j] = -1 (mismatch score) 
 */
int sim_mat[2][2] = {{2, -1}, {-1, 2}};
int gap_score = -2; // gap penalty

// Returns the position of the character in sim_mat
int get_position(char a) {
  int pos;

  if (a == '0') {
    pos = 0;
  } else if (a == '1') {
    pos = 1;
  }
  
  return pos;
}

// Reverse the string s.
string reverse(string s) {
  string reversed;

  for (int i=(int)s.length()-1; i>=0; i--) {
    reversed = reversed + s[i];
  }

  return reversed;
}

// Retutns the number of 1's in s.
int count_ones(string s) {
  int count = 0;
  for (int i=0; i<(int)s.size(); i++) {
    if (s[i] == '1') {
      count++;
    }
  }
  return count;
}

// Returns the number of 1's in y blocks of the first substring.
int count_ones_y(string blocks[], vector<int> y, int pos) {
	int count = 0;
	for (int i=0; i<=pos; i++) {
		count += count_ones(blocks[y[i]]);
	}
	return count;
}

// Returns the number of 1's in x blocks of the second substring.
int count_ones_x(string blocks[], vector<int> x, int pos) {
	int count = 0;
	for (int i=pos; i<(int)x.size(); i++) {
		count += count_ones(blocks[x[i]]);
	}
	return count;
}

// Divide the b_i block into two sets - x-blocks and y-blocks.
void blocks_labeling(string blocks[], int bi[], int bi_size, vector<int> &x,  vector<int> &y) {
  int count_even = 0;
  int count_odd = 0;
  
  for (int i=0; i<bi_size; i++) {
    if (i%2 == 0) {
      count_even += count_ones(blocks[bi[i]]);
    } else {
     count_odd += count_ones(blocks[bi[i]]);
    }
  }
  if (count_even < count_odd) {
    for (int i=0; i<bi_size; i++) {
      if (i%2 == 0) {
				x.push_back(bi[i]);
      } else {
        y.push_back(bi[i]);
      }
    }
  } else {
    for (int i=0; i<bi_size; i++) {
      if (i%2 == 0) {
        y.push_back(bi[i]);
      } else {
        x.push_back(bi[i]);
      }
    }
  }
  
  return;
}
// Returns the sum of the length os sequence of blocks.
int get_length(string blocks[], int pos) {
  int length = 0;
  for (int i=0; i<=pos; i++) {
    length += blocks[i].length();
  }
return length;
}

// Finds the folding point of the protein.
int folding_point(string str) {
  // An array holding the positions of 1's in str.
  int hidrophobic_arr[MAX_SIZE];
  int pos = 0;
  for (int i=0; i < (int)str.length() ;i++) {
    if (str[i] == '1') {
      hidrophobic_arr[pos] = i;
      pos++;
    }
  }
	/* We split the string into blocks. Each block b_i have the form b_i = 1 ot b_i = 1z_01...z_h1, where z_k, 0<=k<=h, is a string 
	 * consisting of only 0 of odd length and h>=1. Then we decompose the string into z_0b_1z_1...b_hz_h, where z_l, 1<=l<=h, is a 
	 * string consisting of only 0 of even length anz z_0 may consist of even or odd length. 
	 */
  string blocks[MAX_SIZE];
  int bi[MAX_SIZE]; // holds the positions of blocks of type b_i
  int bi_iter = 0;
  int j = 0;
  int iter = 0; // the position of the next block in blocks
  // Select z_0 if the string starts with 0.
  if (hidrophobic_arr[0] > 0) {
    blocks[iter] = str.substr(0, hidrophobic_arr[0]);
    iter++;
  }

  bool one = true; // flag showing if we have a block of the type b_i = 1
  while (j+1 < pos) {
    if ((hidrophobic_arr[j+1] - hidrophobic_arr[j]) % 2 == 0) {
      int start = hidrophobic_arr[j];
			// Find the b_i blocks
      while (j+1 < pos && ((hidrophobic_arr[j+1] - hidrophobic_arr[j]) % 2 == 0)) {
        j++;
      }
      int end = hidrophobic_arr[j];
      blocks[iter] = str.substr(start, end-start+1);
      bi[bi_iter] = iter;
      bi_iter++;
      iter++;
      one = false;
    } else {
      if (one) {
        blocks[iter] = "1";
        bi[bi_iter] = iter;
        bi_iter++;
        iter++;
      }
      // Here we have only blocks of 0's - z_i blocks
      if (j < (int)str.length()-1) {
        blocks[iter] = str.substr(hidrophobic_arr[j]+1, hidrophobic_arr[j+1] - hidrophobic_arr[j] - 1); // blokut moje da e prazen, t.e. niama nuli v nego
        iter++;
        j++;
        one = true;
      }
    }
  }

  if (one) {
    blocks[iter] = "1";
    bi[bi_iter] = iter;
    bi_iter++;
    iter++;
  }

  // Select the last block of 0's (if some)
  if (j < (int)str.length()-1) {
    blocks[iter] = str.substr(hidrophobic_arr[j]+1, str.length() - j);
  }

	// Divide the blocks into two cateogies - x-blocks and y-blocks.
  vector<int> x;
	vector<int> y;
  blocks_labeling(blocks, bi, bi_iter, x, y);
	bool y_is_first;// flag
	if (x[0] < y[0]) {
		y_is_first = 0;
	} else {
		y_is_first = 1;
	}
	vector<int>::iterator block_pos_x;
  vector<int>::iterator block_pos_y;
  vector<int> min_ones;
  vector<int> lengths;
  int left_length = 0;

	for (int i=0; i<bi_iter; i++) {
		if ((block_pos_y = find(y.begin(), y.end(), bi[i])) != y.end()) {
			int y_pos = block_pos_y - y.begin();
      int count_y = count_ones_y(blocks, y, y_pos);
			int x_pos = (y_is_first) ? y_pos : (y_pos+1);
      if (x_pos < x.size() && x_pos >= 0) {
        int count_x = count_ones_x(blocks, x, x_pos);
        min_ones.push_back(min(count_x, count_y));
        lengths.push_back(get_length(blocks, bi[i]));
      }
		} else if ((block_pos_x = find(x.begin(), x.end(), bi[i])) != x.end()) {
			int x_pos = block_pos_x - x.begin();
			if (x_pos !=0 ) {
        int count_x = count_ones_x(blocks, x, x_pos+1);
				int y_pos = (y_is_first) ? x_pos : x_pos-1;
        if (y_pos < y.size() && y_pos >= 0) {
          int count_y = count_ones_y(blocks, y, y_pos);
          min_ones.push_back(min(count_x, count_y));
          lengths.push_back(get_length(blocks, bi[i]));
        }
      }
		}
	}

	int max_pos = max_element(min_ones.begin(),min_ones.end()) - min_ones.begin();
	return lengths[max_pos];
}

void alignment(string s, string t, string &s_aln, string &t_aln) {
	// Sequence Alignment - if we have a gap in one of the two substring, we add two gaps so that the number of gaps is always with even length
  int n = (int)s.length();
  int m = (int)t.length();
  int D[MAX_SIZE][MAX_SIZE];
  memset(D, 0, sizeof(D));
  D[1][1] = 0;

  for (int i=1; i<=n; i++) {
   D[1][i+1] = gap_score*i;
  }

  for (int j=1; j<=m; j++) {
   D[j+1][1] = gap_score*j;
  }

  for (int i = 2; i <= m+1; i++) {
    for (int j = 2; j<= n+1; j++) {
      int match = D[i-1][j-1] + sim_mat[get_position(s[j-2])][get_position(t[i-2])];
      int gaps = D[i][j-1] + gap_score;
      int gapt = D[i-1][j] + gap_score;
      int max = match;
      max = (gaps > max) ? gaps : max;
      max = (gapt > max) ? gapt : max;
      D[i][j] = max;
    }
  }

  int i = m+1;
  int j = n+1;
  while( i > 1 && j > 1) {
    if ( (D[i][j] - sim_mat[get_position(s[j-2])][get_position(t[i-2])] == D[i-1][j-1]) && t[i-2] == s[j-2]) {
      t_aln = t[i-2] + t_aln;
      s_aln = s[j-2] + s_aln;
      i = i-1; j = j-1;
    } else if ( (D[i][j] - sim_mat[get_position(s[j-2])][get_position(t[i-2])] == D[i-1][j-1]) && t[i-2] != s[j-2]) {
      if (s[j-2] == '0') {
        s_aln = s[j-2] + s_aln;
				t_aln = "-" + t_aln;
				if (j > 2) {
					s_aln = s[j-3] + s_aln;
					t_aln = "-" + t_aln;
					j--;
				}
				j--;
			} else if (t[i-2] == '0') {
        t_aln = t[i-2] + t_aln;
				s_aln = "-" + s_aln;
				if (i > 2) {
					t_aln = t[i-3] + t_aln;
					s_aln = "-" + s_aln;
					i--;
				}
				i--;
			}
		} else if ( (D[i][j] - gap_score) == D[i-1][j]) {
      s_aln = "-" + s_aln;
      t_aln = t[i-2] + t_aln;
			if (i > 2) {
				t_aln = t[i-3] + t_aln;
				s_aln = "-" +s_aln;
				i--;
			}
			i--;
		} else if ((D[i][j] - gap_score) == D[i][j-1]) {
      s_aln = s[j-2] + s_aln;
			t_aln = "-" + t_aln;
			if (j > 2) {
				s_aln = s[j-3] + s_aln;
				t_aln = "-" + t_aln;
				j--;
			}
			j--;
    }
  }

 // Fill the rest of the shorter substring with gaps. 
  if (j > 1) {
    while (j > 1) {
      s_aln = s[j-2] + s_aln;
      t_aln = '-' + t_aln;
      j = j-1;
    }
  } else if (i > 1) {
    while (i > 1) {
      s_aln = '-' + s_aln;
      t_aln = t[i-2] + t_aln;
      i = i-1;
    }
  }

	// Fix alignments where the block of 0 is separated by gaps, or we have 1 surrounded by gaps.
  int ind = 0;
  while (ind < (int) s_aln.length()) {
    if (s_aln[ind] == t_aln[ind]) {
      ind++;
    }
    if (s_aln[ind] == '-') {
      int cnt = 0;
      while (s_aln[ind] == '-') {
        cnt++;
        ind++;
      }
			if ((s_aln[ind] == '0' && t_aln[ind] == '0') || (s_aln[ind] == '1' && t_aln[ind] == '1' && s_aln[ind+1] == '-')) {
        int new_ind = ind - cnt;
				if (new_ind > 0) { 
					while (s_aln[new_ind] == '-') {
						new_ind--;
					}
					new_ind++;
				}
        if (t_aln[new_ind] == t_aln[ind]) {
					if (t_aln[ind] == '0' && t_aln[new_ind+1] != '1') {
						s_aln[new_ind] = t_aln[ind];
						s_aln[ind] = '-';
					}
        }
			}
    } else if (t_aln[ind] == '-') {
      int cnt = 0;
      while (t_aln[ind] == '-') {
        cnt++;
        ind++;
      }
			if ((t_aln[ind] == '0' && s_aln[ind] == '0') || (t_aln[ind] == '1' && s_aln[ind] == '1' && t_aln[ind+1] == '-')) {
        int new_ind = ind - cnt;
				if (new_ind > 0) { 
					while (t_aln[new_ind] == '-') {
						new_ind--;
					}
					new_ind++;
				}
        if (s_aln[new_ind] == s_aln[ind]) {
					if (s_aln[ind] == '0' && t_aln[new_ind+1] != '1') {
						t_aln[new_ind] = s_aln[ind];
						t_aln[ind] = '-';
					}
        }
			}
    }
  }
}

void print_protein (string s, string t) {
  char lattice[LAT_SIZE][LAT_SIZE];
  memset(lattice, ' ', sizeof(lattice));
  int pos = MAX_SIZE/2;
  int size = (int)s.length();
  int j = 0;
  int starting_gaps = 0;
	/* If one of the aligned substring starts with gaps, we just print only the one without gaps */
	if (s[starting_gaps] == '-') {
		while (s[starting_gaps] == '-') {
			starting_gaps++;
		}
		for (int i=0; i<starting_gaps; i++) {
			lattice[pos][j] = t[i];
			lattice[pos][j+1] = '-';
			j+=2;
		}
	} else	if (t[starting_gaps] == '-') {
		while (t[starting_gaps] == '-') {
			starting_gaps++;
		}
		for (int i=0; i<starting_gaps; i++) {
			lattice[pos+2][j] = s[i];
			lattice[pos+2][j+1] = '-';
			j+=2;
		}
	}

	/* If we have gaps at the and of one of the aligned substrings, we fold the protein at the middle of the gaps. 
	 * We have constructed the alignment so that all the gaps have an even length. What's more these gaps corespond
	 * to a z type of block. */
  int s_ending_gaps = 0 ;
  int t_ending_gaps = 0;
  if (s[s.length()-1] == '-') {
		int gap_pos = (int)s.length() - 1;
		while (s[gap_pos] == '-') {
			s_ending_gaps++;
			gap_pos--;
		}
  }

  if (t[s.length()-1] == '-') {
    int gap_pos = (int)t.length()-1;
    while (t[gap_pos] == '-') {
      t_ending_gaps++;
      gap_pos--;
    }
  }

  int length = (int)s.length() - max(s_ending_gaps, t_ending_gaps);
  for (int i = starting_gaps; i < (int)length; i++) {
		if (s[i] != '-' && t[i] != '-') {
      lattice[pos][j] = t[i];
      lattice[pos+2][j] = s[i];
      if (s[i+1] && t[i+1] && s[i+1] != '-' && t[i+1] != '-') {
        lattice[pos][j+1] = '-';
        lattice[pos+2][j+1] = '-';
        j++;
      }
      j++;
		} else if (s[i] == '-') {
			if (i == starting_gaps) {
				j--;
			}
      int cnt = 0;
      while (s[i] == '-') {
				i++;
				cnt++;
      }
      for (int a=1; a<=cnt/2; a++) {
        lattice[pos-2*a][j-1] = t[i-cnt-1+a];
        lattice[pos-2*a][j+1] = t[i-a]; 
				if (i-cnt != starting_gaps) {
					lattice[pos-2*a+1][j-1] = '|';
				}
				lattice[pos-2*a+1][j+1] = '|';
      }
      lattice[pos-cnt][j] = '-';
			lattice[pos+2][j] = '-';
      j+=1;
      i--;
    } else if(t[i] == '-') {
			if (i == starting_gaps) {
				j--;
			}
      int cnt = 0;
      while (t[i] == '-') {
        i++;
        cnt++;
      }
      for (int a=1; a<=cnt/2; a++) {
        lattice[pos + 2*a + 2][j-1] = s[i-cnt-1+a];
        lattice[pos + 2*a + 2 ][j+1] = s[i-a];
				if (i-cnt != starting_gaps) {
					lattice[pos + 2*a + 1][j-1] = '|';
				}
				lattice[pos + 2*a + 1][j+1] = '|';
      }
      lattice[pos+2+cnt][j] = '-';
      lattice[pos][j] = '-';
      j += 1;
      i--;
    }
  }

  if (t_ending_gaps > 0 && t_ending_gaps%2 == 0) {
    lattice[pos][j] = '-';
    lattice[pos+2][j] = '-';
    j++;
    int t_last = (int)t.length()-1;
    int t_first = t_last - t_ending_gaps + 1;
    int half_gaps = t_ending_gaps/2;
    for (int i=0; i < half_gaps; i++) {
      lattice[pos][j] = s[t_last-i];
      lattice[pos+2][j] = s[t_first+i];
      if (i != half_gaps-1) {
        lattice[pos][j+1] = '-';
        lattice[pos+2][j+1] = '-';
        j++;
      }
      j++;
    }
  }

  if (s_ending_gaps > 0 && s_ending_gaps%2 == 0) {
    lattice[pos][j] = '-';
    lattice[pos+2][j] = '-';
    j++;
    int s_last = (int)s.length()-1;
    int s_first = s_last - s_ending_gaps + 1;
    int half_gaps = s_ending_gaps/2;
    for (int i=0; i < half_gaps; i++) {
      lattice[pos][j] = t[s_last-i];
      lattice[pos+2][j] = t[s_first+i];
      if (i != half_gaps-1) {
        lattice[pos][j+1] = '-';
        lattice[pos+2][j+1] = '-';
        j++;
      }
    }
    j++;
  }

  lattice[pos+1][j-1] = '|';
  for (int i=0 ; i<MAX_SIZE; i++) {
    for (int j=0; j<MAX_SIZE; j++) {
      cout<<lattice[i][j];
    }
    cout<<endl;
  }
}

int main(){
	//string str = "0100101001110101000010"; //ok
  //string str = "00001100100100100100101000000100"; //ok
  //string str = "0100001010011000011010"; //print issues
  //string str = "010010000001000110101000101100000010000100000000010"; //print issues
  //string str = "11100010100001000000100101010000010"; //ok
	//string str = "01001010000010101001000100101010010"; //print issues
	//string str = "00100011100001001010101010101"; //ok
	//string str = "0010010100101000100100100000101001011"; - no
	string str = "1000010101010101011011101010110";
	int folding_point_pos = folding_point(str);
  string t = str.substr(0, (size_t)folding_point_pos);
  string right = str.substr((size_t)folding_point_pos, str.length());
  string s = reverse(right);
  string t_aln;
  string s_aln;
	cout<<"str = "<<str<<endl;
  cout<<"t = "<<t<<endl;
  cout<<"s = "<<s<<endl;
  alignment(s, t, s_aln, t_aln);
	cout<<"t aligned is "<<t_aln<<endl;
  cout<<"s aligned is "<<s_aln<<endl;
  print_protein (s_aln, t_aln);

  return 0;
}