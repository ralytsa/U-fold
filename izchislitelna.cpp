#include <iostream>
#include <sstream>
#include <string>
#include <hash_map>
#include <cstring>

using namespace std;

#define MAX_SIZE 50
#define LAT_SIZE 100

int sim_mat[2][2] = {{2, -1}, {-1, 2}};
int gap_score = -2;

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

// Retutns the number of 1's int s.
int count_ones(string s) {
  int count = 0;

  for (int i=0; i<(int)s.size(); i++) {
    if (s[i] == '1') {
      count++;
    }
  }

  return count;
}

// Returns the difference 1's in the two sets.
// The number of 1's in the two sets must be almost equal.
int folding_point_set_difference(string blocks[], int size, int pos) {
  int left_count_ones = 0;
  int right_count_ones = 0;

  for (int i=0; i<=pos; i++) {
    left_count_ones += count_ones(blocks[i]);
  }
  for(int i=pos+1; i<size; i++) {
    right_count_ones += count_ones(blocks[i]);
  }

  return right_count_ones - left_count_ones;
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

  string blocks[MAX_SIZE];
  int j = 0;
  int iter = 0; // the position of the next block in blocks

  if (hidrophobic_arr[0] > 0) {
    blocks[iter] = str.substr(0, hidrophobic_arr[0]);
    iter++;
  }

  bool one = true; // flag showing if we have a block of the type b_i = 1
  while (j+1 < pos) {
    if ((hidrophobic_arr[j+1] - hidrophobic_arr[j]) % 2 == 0) {
      int start = hidrophobic_arr[j];
      while (j+1 < pos && ((hidrophobic_arr[j+1] - hidrophobic_arr[j]) % 2 == 0)) {
        j++;
      }
      int end = hidrophobic_arr[j];
      blocks[iter] = str.substr(start, end-start+1);
      iter++;
      one = false;
    } else {
      if (one) {
        blocks[iter] = "1";
        iter++;
      }
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
    iter++;
  }

  if (j < (int)str.length()-1) {
    blocks[iter] = str.substr(hidrophobic_arr[j]+1, str.length() - j);
  }

  int left_length = 0;
  for (int i = 0; i <= iter; i++) {
    if ((int)blocks[i].length() > 0) {
      int distance = folding_point_set_difference(blocks, iter+1, i);
      if (distance > 0) {
        left_length += (int)blocks[i].length();
      } else {
        int old_distance = folding_point_set_difference(blocks, iter+1, i-1);
        if (abs(distance) < abs(old_distance)) {
          left_length += (int)blocks[i].length();
        } else {
          return left_length;
        }
      }
    }
  }

  return 0;
}

void print_protein (string s, string t) {
  char lattice[LAT_SIZE][LAT_SIZE];
  memset(lattice, ' ', sizeof(lattice));
  int pos = MAX_SIZE/2;
  int size = (int)s.length();
  int j = 0;
  int starting_zeros = 0;
  if ((s[starting_zeros] == '0' || s[starting_zeros] == '-') && (t[starting_zeros] == '0' || t[starting_zeros] == '-')) {
    while ((s[starting_zeros] == '0' || s[starting_zeros] == '-') && (t[starting_zeros] == '0' || t[starting_zeros] == '-')) {
      if (t[starting_zeros] == '0') {
        lattice[pos][j] = t[starting_zeros];
        lattice[pos][j+1] = '-';
      }
      if (s[starting_zeros] == '0') {
        lattice[pos+2][j] = s[starting_zeros];
        lattice[pos+2][j+1] = '-';
      }
      j+=2;
      starting_zeros++;
    }
    if (s[starting_zeros] == '-' || t[starting_zeros] == '-') {
    j--;
    lattice[pos][j] = ' ';
    lattice[pos+2][j] = ' ';
    }
  }
  int s_ending_gaps = 0 ;
  int t_ending_gaps = 0;
  if (s[s.length()-1] == '-') {
  int gap_pos = s.length() - 1;
  while (s[gap_pos] == '-') {
    s_ending_gaps++;
    gap_pos--;
  }
  }

  if (t[s.length()-1] == '-') {
    int gap_pos = t.length()-1;
    while (t[gap_pos] == '-') {
      t_ending_gaps++;
      gap_pos--;
    }
  }
  int length = (int)s.length() - max(s_ending_gaps, t_ending_gaps);
  for (int i = starting_zeros; i < (int)length; i++) {
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
      int cnt = 0;
      while (s[i] == '-') {
      i++;
      cnt++;
      }
      for (int a=1; a<=cnt/2; a++) {
        lattice[pos-2*a][j-1] = t[i-cnt/2-a];
        lattice[pos-2*a][j+1] = t[i-a];
        lattice[pos-2*a+1][j-1] = '|';
        lattice[pos-2*a+1][j+1] = '|';
      }
      lattice[pos-cnt][j] = '-';
      lattice[pos+2][j] = '-';
      j+=1;
      i--;
    } else if(t[i] == '-') {
      int cnt = 0;
      while (t[i] == '-') {
        i++;
        cnt++;
      }
      for (int a=1; a<=cnt/2; a++) {
        lattice[pos + 2*a + 2][j-1] = s[i-cnt/2-a];
        lattice[pos + 2*a + 2 ][j+1] = s[i-a];
        lattice[pos + 2*a + 1][j-1] = '|';
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


string reverse(string s) {
  string reversed;

  for (int i=(int)s.length()-1; i>=0; i--) {
    reversed = reversed + s[i];
  }

  return reversed;
}

void alignment(string s_org, string t_org, string &s_aln, string &t_aln) {
  int s_starting_zeros = 0;
  while (s_org[s_starting_zeros] == '0') {
   s_starting_zeros++;
  }
  string s = s_org.substr((size_t)s_starting_zeros, s_org.length());
  int t_starting_zeros = 0;
  while (t_org[t_starting_zeros] == '0') {
   t_starting_zeros++;
  }
  string t = t_org.substr((size_t)t_starting_zeros, t_org.length());

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
        s_aln = s[j-3] + s_aln;
        t_aln = "--" + t_aln;
        j = j-2;
      } else if (t[i-2] == '0') {
        t_aln = t[i-2] + t_aln;
        t_aln = t[i-3] + t_aln;
        s_aln = "--" + s_aln;
        i = i-2;
      }
    } else if ( (D[i][j] - gap_score) == D[i-1][j]) {
      s_aln = "--" + s_aln;
      t_aln = t[i-2] + t_aln;
      t_aln = t[i-3] + t_aln;
      i = i-2;
    } else if ((D[i][j] - gap_score) == D[i][j-1]) {
      s_aln = s[j-2] + s_aln;
      s_aln = s[j-3] + s_aln;
      t_aln = "--" + t_aln;
      j = j-2;
    }
  }

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
  cout<<t_aln<<endl;
  cout<<s_aln<<endl;
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
      if (s_aln[ind] == '0' && t_aln[ind] == '0' ) {
        int new_ind = ind - cnt;
        while (s_aln[new_ind] == '-') {
          new_ind--;
        }
        if (t_aln[new_ind+1] == '0') {
          s_aln[new_ind+1] = '0';
          s_aln[ind] = '-';
        }
      }
    }
    if (t_aln[ind] == '-') {
      int cnt = 0;
      while (t_aln[ind] == '-') {
        cnt++;
        ind++;
      }
      if (t_aln[ind] == '0' && s_aln[ind] == '0' ) {
        int new_ind = ind - cnt;
        while (t_aln[new_ind] == '-') {
          new_ind--;
        }
        if (s_aln[new_ind+1] == '0') {
          t_aln[new_ind+1] = '0';
          t_aln[ind] = '-';
        }
      }
    }
  }

  if (s_starting_zeros >= t_starting_zeros) {
    for (int i=0; i<t_starting_zeros; i++) {
      s_aln = "0" + s_aln;
      t_aln = "0" + t_aln;
    }
    int diff = s_starting_zeros - t_starting_zeros;
    for (int i=0; i<diff; i++) {
      s_aln = "0" + s_aln;
      t_aln = "-" + t_aln;
    }
  } else {
    for (int i=0; i<s_starting_zeros; i++) {
      s_aln = "0" + s_aln;
      t_aln = "0" + t_aln;
    }
    int diff = t_starting_zeros - s_starting_zeros;
    for (int i=0; i<diff; i++) {
      t_aln = "0" + t_aln;
      s_aln = "-" + s_aln;
    }
  }
}


int main(){
  //string str = "0100101001110101000010";
  //string str = "00001100100100100100101000000100";
  //string str = "0100001010011000011010";
  //string str = "010010000001000110101000101100000010000100000000010"; 
  string str = "11100010100001000000100101010000010";
  int folding_point_pos = folding_point(str);
  string t = str.substr(0, (size_t)folding_point_pos);
  string right = str.substr((size_t)folding_point_pos, str.length());
  string s = reverse(right);
  string t_aln;
  string s_aln;
  cout<<t<<endl;
  cout<<s<<endl;
  alignment(s, t, s_aln, t_aln);
  cout<<t_aln<<endl;
  cout<<s_aln<<endl;
  print_protein (s_aln, t_aln);
}