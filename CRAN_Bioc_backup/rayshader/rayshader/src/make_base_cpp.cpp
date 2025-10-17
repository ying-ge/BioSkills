#include <Rcpp.h>
using namespace Rcpp;

template <int RTYPE>
Matrix<RTYPE> vec2matrix(Vector<RTYPE> x, int nrow, int ncol) {
  return Matrix<RTYPE>(nrow, ncol, x.begin());
}

// [[Rcpp::export]]
List make_surface_cpp(NumericMatrix& heightmap,
                      LogicalMatrix& na_matrix,
                      NumericMatrix& normalsx,
                      NumericMatrix& normalsy,
                      NumericMatrix& normalsz,
                      double basedepth) {
  std::vector<NumericMatrix> vertices;
  std::vector<NumericMatrix> normals;
  std::vector<NumericMatrix> texcoords;
  
  int rows = heightmap.nrow();
  int cols = heightmap.ncol();
  for(int j = 0; j < rows-1; j++) {
    for(int i = 0; i < cols-1; i++) {
      if(!na_matrix(j,i) && !na_matrix(j, i + 1) && !na_matrix(j+1, i) && !na_matrix(j+1,i+1)) {
        vertices.push_back(vec2matrix(NumericVector::create(j,j+1,j, 
                                                            heightmap(j,i),heightmap(j+1,i),heightmap(j,i+1), 
                                                            i,i,i+1),3,3));
        vertices.push_back(vec2matrix(NumericVector::create(j+1,j+1,j, 
                                                            heightmap(j+1,i),heightmap(j+1,i+1),heightmap(j,i+1), 
                                                            i,i+1,i+1),3,3));
        normals.push_back(vec2matrix(NumericVector::create(normalsx(j,i),normalsx(j+1,i),normalsx(j,i+1), 
                                                           normalsy(j,i),normalsy(j+1,i),normalsy(j,i+1), 
                                                           normalsz(j,i),normalsz(j+1,i),normalsz(j,i+1)),3,3));
        normals.push_back(vec2matrix(NumericVector::create(normalsx(j+1,i),normalsx(j+1,i+1),normalsx(j,i+1), 
                                                           normalsy(j+1,i),normalsy(j+1,i+1),normalsy(j,i+1), 
                                                           normalsz(j+1,i),normalsz(j+1,i+1),normalsz(j,i+1)),3,3));
        texcoords.push_back(vec2matrix(NumericVector::create((float)j/(float)rows, (float)(j+1)/(float)rows,(float)j/(float)rows,
                                                             (float)i/(float)cols,(float)i/(float)cols,(float)(i+1)/(float)cols),3,2));
        texcoords.push_back(vec2matrix(NumericVector::create((float)(j+1)/(float)rows,(float)(j+1)/(float)rows,(float)j/(float)rows,
                                                               (float)i/(float)cols,(float)(i+1)/(float)cols,(float)(i+1)/(float)cols),3,2));
      } 
    }
  }
  List vectorlist = wrap(vertices);
  List normallist = wrap(normals);
  List texcoordlist = wrap(texcoords);
  
  return(List::create(_["vertices"] = vectorlist, _["normals"] = normallist, _["texcoords"] = texcoordlist));
}

// [[Rcpp::export]]
List make_base_cpp(NumericMatrix& heightmap,
                    LogicalMatrix& na_matrix,
                    double basedepth) {
  std::vector<NumericMatrix> vertices;
  // std::vector<NumericMatrix> normals;
  std::vector<bool> horizontal;
  std::vector<double> height;
  
  int rows = heightmap.nrow();
  int cols = heightmap.ncol();
  for(int j = 0; j < rows-1; j++) {
    for(int i = 0; i < cols-1; i++) {
      if(j == 0) {
        if(!std::isnan(heightmap(0+j,i)) && !std::isnan(heightmap(0+j,i+1))) {
          if(!na_matrix(j,i) && !na_matrix(j, i + 1)) {
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(0+j,i),basedepth,basedepth, -i-1,-i-2,-i-1),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(0+j,i),heightmap(0+j,i+1),basedepth, -i-1,-i-2,-i-2),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(-1,-1,-1, 0,0,0, 0,0,0),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(-1,-1,-1, 0,0,0, 0,0,0),3,3));
            horizontal.push_back(true);
            horizontal.push_back(true);
            height.push_back( heightmap(0+j,i) - basedepth);
            height.push_back( heightmap(0+j,i+1) - basedepth);

          }
        }
      } else {
        if(!std::isnan(heightmap(0+j,i)) && !std::isnan(heightmap(0+j,i+1))) {
          if((!na_matrix(j,i) && na_matrix(j - 1, i) && !na_matrix(j, i+1)) || (!na_matrix(j,i+1) && na_matrix(j - 1, i+1))) {
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(0+j,i),basedepth,basedepth, -i-1,-i-2,-i-1),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(0+j,i),heightmap(0+j,i+1),basedepth, -i-1,-i-2,-i-2),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(-1,-1,-1, 0,0,0, 0,0,0),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(-1,-1,-1, 0,0,0, 0,0,0),3,3));
            horizontal.push_back(true);
            horizontal.push_back(true);
            height.push_back( heightmap(0+j,i)  - basedepth);
            height.push_back( heightmap(0+j,i+1)  - basedepth);
          }
        }
      }
    }
  }
  
  for(int j = 0; j < rows; j++) {
    for(int i = 0; i < cols-1; i++) {
      if(j == rows - 1) {
        if(!std::isnan(heightmap(j,i)) && !std::isnan(heightmap(j,i+1))) {
          if(!na_matrix(j,i+1)) {
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(j,i),basedepth,basedepth, -i-1, -i-1, -i-2),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(j,i),basedepth,heightmap(j,i+1), -i-1,-i-2,-i-2),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(1,1,1, 0,0,0, 0,0,0),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(1,1,1, 0,0,0, 0,0,0),3,3));
            horizontal.push_back(true);
            horizontal.push_back(true);
            height.push_back( heightmap(j,i)  - basedepth);
            height.push_back( heightmap(j,i+1)  - basedepth);
          }
        }
      } else {
        if(!std::isnan(heightmap(j,i)) && !std::isnan(heightmap(j,i+1))) {
          if((!na_matrix(j,i) && na_matrix(j+1, i) && !na_matrix(j, i+1)) || (!na_matrix(j,i+1) && na_matrix(j+1, i+1))) {
            vertices.push_back(vec2matrix(NumericVector::create(j+1,j+1,j+1, heightmap(j,i),basedepth,basedepth, -i-1, -i-1, -i-2),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(j+1,j+1,j+1, heightmap(j,i),basedepth,heightmap(j,i+1), -i-1,-i-2,-i-2),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(1,1,1, 0,0,0, 0,0,0),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(1,1,1, 0,0,0, 0,0,0),3,3));
            horizontal.push_back(true);
            horizontal.push_back(true);
            height.push_back( heightmap(0+j,i)  - basedepth);
            height.push_back( heightmap(0+j,i+1)  - basedepth);
          }
        }
      }
    }
  }
  for(int j = 0; j < cols-1; j++) {
    for(int i = 0; i < rows-1; i++) {
      if(j == 0) {
        if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
          if(!na_matrix(i+1,j)) {
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+1,i+2, heightmap(i,0+j),basedepth,basedepth,  -1+j,-1+j,-1+j),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+2, heightmap(i,0+j),basedepth,heightmap(i+1,0+j), -1+j,-1+j,-1+j),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,-1,-1,-1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,-1,-1,-1),3,3));
            horizontal.push_back(false);
            horizontal.push_back(false);
            height.push_back( heightmap(i,j) - basedepth);
            height.push_back( heightmap(i+1,j) - basedepth);
          }
        }
      } else {
        if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
          if((!na_matrix(i,j) && na_matrix(i,j - 1) && !na_matrix(i+1,j)) || (!na_matrix(i,j) && na_matrix(i + 1, j-1) && !na_matrix(i+1,j))) {
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+1,i+2, heightmap(i,0+j),basedepth,basedepth,  -1-j,-1-j,-1-j),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+2, heightmap(i,0+j),basedepth,heightmap(i+1,0+j), -1-j,-1-j,-1-j),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,-1,-1,-1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,-1,-1,-1),3,3));
            horizontal.push_back(false);
            horizontal.push_back(false);
            height.push_back( heightmap(i,j) - basedepth);
            height.push_back( heightmap(i+1,j) - basedepth);
          }
        }
      }
    }
  }
  for(int j = 0; j < cols; j++) {
    for(int i = 0; i < rows-1; i++) {
      if(j == cols - 1) {
        if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
          if(!na_matrix(i,j) && !na_matrix(i+1,j)) {
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+1, heightmap(i,j),basedepth,basedepth,-j-1,-j-1,-j-1),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+2, heightmap(i,j),heightmap(i+1,j),basedepth, -j-1,-j-1,-j-1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,1,1,1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,1,1,1),3,3));
            horizontal.push_back(false);
            horizontal.push_back(false);
            height.push_back( heightmap(i,j) - basedepth);
            height.push_back( heightmap(i+1,j) - basedepth);
          }
        }
      } else {
        if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
          if((!na_matrix(i,j) && na_matrix(i,j + 1) && !na_matrix(i+1,j)) || (!na_matrix(i,j) && na_matrix(i+1, j+1) && !na_matrix(i+1,j))) {
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+1, heightmap(i,j),basedepth,basedepth,-j-1,-j-1,-j-1),3,3));
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+2,i+2, heightmap(i,j),heightmap(i+1,j),basedepth, -j-1,-j-1,-j-1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,1,1,1),3,3));
            // normals.push_back(vec2matrix(NumericVector::create(0,0,0, 0,0,0,1,1,1),3,3));
            horizontal.push_back(false);
            horizontal.push_back(false);
            height.push_back( heightmap(i,j) - basedepth);
            height.push_back( heightmap(i+1,j) - basedepth);
          }
        }
      }
    }
  }
  List vectorlist = wrap(vertices);
  // List normallist = wrap(normals);
  LogicalVector directionlist = wrap(horizontal);
  NumericVector heights = wrap(height);
  
  return(List::create(_["vertices"] = vectorlist, //_["normals"] = normallist, 
                      _["is_horizontal"] = directionlist, _["edge_heights"] = heights));
}

// [[Rcpp::export]]
List make_water_cpp(NumericMatrix& heightmap,
                    LogicalMatrix& na_matrix,
                    double waterheight) {
  int rows = heightmap.nrow();
  int cols = heightmap.ncol();
  std::vector<NumericMatrix> vertices;
  double endcoord, begincoord, heighttemp;
  int offset = 1;
  double adjust;
  for(int j = 0; j < rows - 1; j++) {
    offset = 0;
    if(j != 0) {
      offset = 1;
    }
    for(int i = 0; i < cols - 1; i++) {
      if(!std::isnan(heightmap(j,i)) && !std::isnan(heightmap(j,i+1))) {
        if(((na_matrix(j-offset,i) && !na_matrix(j,i) && !na_matrix(j,i+1)) || (na_matrix(j-offset,i+1) && !na_matrix(j,i+1))) || j == 0) {
          if(heightmap(0+j,i) < waterheight && heightmap(0+j,i+1) < waterheight) {
            adjust = (waterheight - heightmap(0+j,i))/(heightmap(0+j,i+1)-heightmap(0+j,i));
            if(heightmap(0+j,i+1) > waterheight && fabs(adjust) < 1) {
              endcoord = -(double)i - adjust;
            } else {
              endcoord = -i - 1;
            }
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heightmap(0+j,i),waterheight,waterheight, -i,-i,endcoord),3,3));
            if(heightmap(0+j,i) > waterheight && fabs(adjust) < 1) {
              begincoord = -(double)i - adjust;
              heighttemp = waterheight;
            } else {
              begincoord = -i;
              heighttemp = heightmap(0+j,i);
            }
            vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,1+j, heighttemp,waterheight,heightmap(0+j,i+1), begincoord,-i-1,-i-1),3,3));
          }
        }
      }
    }
  }
  for(int j = 0; j < cols - 1; j++) {
    offset = 0;
    if(j != 0) {
      offset = 1;
    }
    for(int i = 0; i < rows-1; i++) {
      if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
        if(((na_matrix(i,j-offset) && !na_matrix(i,j) && !na_matrix(i+1,j)) || (na_matrix(i+1,j-offset) && !na_matrix(i+1,j))) || j == 0) {
          if(heightmap(i,0+j) < waterheight && heightmap(i+1,0+j) < waterheight) {
            adjust = (waterheight - heightmap(i,0+j))/(heightmap(i+1,0+j)-heightmap(i,0+j));
            if(heightmap(i+1,0+j) > waterheight && fabs(adjust) < 1) {
              endcoord = (double)i + 1 + adjust;
            } else {
              endcoord = i+2;
            }
            vertices.push_back(vec2matrix(NumericVector::create(i+1,endcoord,i+1, heightmap(i,0+j),waterheight,waterheight,  -j,-j,-j),3,3));
            if(heightmap(i,0+j) > waterheight && fabs(adjust) < 1) {
              begincoord = (double)i + 1 + adjust;
              heighttemp = waterheight;
            } else {
              begincoord = i+1;
              heighttemp = heightmap(i,0+j);
            }
            vertices.push_back(vec2matrix(NumericVector::create(begincoord,i+2,i+2, heighttemp,heightmap(i+1,0+j),waterheight, -j,-j,-j),3,3));
          }
        }
      }
    }
  }
  for(int j = 0; j < rows; j++) {
    offset = 0;
    if(j != rows - 1) {
      offset = 1;
    }
    for(int i = 0; i < cols - 1; i++) {
      if(!std::isnan(heightmap(j,i)) && !std::isnan(heightmap(j,i+1))) {
        if(((na_matrix(j+offset,i) && !na_matrix(j,i) && !na_matrix(j,i+1)) || (na_matrix(j+offset,i+1) && !na_matrix(j,i+1))) || j == rows-1) {
          if(heightmap(j,i) < waterheight && heightmap(j,i+1) < waterheight) {
            adjust = (waterheight - heightmap(j,i))/(heightmap(j,i+1)-heightmap(j,i));
            if(heightmap(j,i+1) > waterheight && fabs(adjust) < 1) {
              endcoord = -(double)i - adjust;
            } else {
              endcoord = -i - 1;
            }
            vertices.push_back(vec2matrix(NumericVector::create(j+1,j+1,j+1, heightmap(j,i),waterheight,waterheight, -i, endcoord, -i),3,3));
            if(heightmap(j,i) > waterheight && fabs(adjust) < 1) {
              begincoord = -(double)i - adjust;
              heighttemp = waterheight;
            } else {
              begincoord = -i;
              heighttemp = heightmap(j,i);
            }
            vertices.push_back(vec2matrix(NumericVector::create(j+1,j+1,j+1, heighttemp,heightmap(j,i+1),waterheight, begincoord,-i-1,-i-1),3,3));
          }
        }
      }
    }
  }
  for(int j = 0; j < cols; j++) {
    offset = 0;
    if(j != cols - 1) {
      offset = 1;
    }
    for(int i = 0; i < rows-1; i++) {
      if(!std::isnan(heightmap(i,j)) && !std::isnan(heightmap(i+1,j))) {
        if(((na_matrix(i,j+offset) && !na_matrix(i,j) && !na_matrix(i+1,j)) || (na_matrix(i+1,j+offset) && !na_matrix(i+1,j))) || j == cols-1) {
          if(heightmap(i,j) < waterheight && heightmap(i+1,j) < waterheight) {
            adjust = (waterheight - heightmap(i,j))/(heightmap(i+1,j)-heightmap(i,j));
            if(heightmap(i+1,j) > waterheight && fabs(adjust) < 1) {
              endcoord = (double)i + adjust;
            } else {
              endcoord = i+2;
            }
            vertices.push_back(vec2matrix(NumericVector::create(i+1,i+1,endcoord, heightmap(i,j),waterheight,waterheight,-j,-j,-j),3,3));
            if(heightmap(i,j) > waterheight && fabs(adjust) < 1) {
              begincoord = (double)i  + adjust;
              heighttemp = waterheight;
            } else {
              begincoord = i+1;
              heighttemp = heightmap(i,j);
            }
            vertices.push_back(vec2matrix(NumericVector::create(begincoord,i+2,i+2,  heighttemp,waterheight,heightmap(i+1,j), -j,-j,-j),3,3));
          }
        }
      }
    }
  }
  List vectorlist = wrap(vertices);
  return(vectorlist);
}

// [[Rcpp::export]]
List make_waterlines_cpp(NumericMatrix& heightmap,
                        LogicalMatrix& na_matrix,
                         double waterdepth) {
  std::vector<NumericMatrix> vertices;
  int rows = heightmap.nrow();
  int cols = heightmap.ncol();
  int offset, offset2 = 0;
  int offsetside, offsetside2 = 0;
  bool drawing = false;
  double startcoord, endcoord = 1;
  double adjust;
  
  for(int j = 0; j < rows; j++) {
    drawing = false;
    if(j != 0) {
      offset = 1;
    } else {
      offset = 0;
    }
    if(j != rows-1) {
      offset2 = 1;
    } else {
      offset2 = 0;
    }
    for(int i = 0; i < cols; i++) {
      if(i != 0) {
        offsetside = 1;
      } else {
        offsetside = 0;
      }
      if(i != cols-1) {
        offsetside2 = 1;
      } else {
        offsetside2 = 0;
      }
      //Edges
      if(drawing && (j == 0 || j == rows - 1)) {
        if((heightmap(j,i) > waterdepth || i == cols-1 ) || na_matrix(j,i)) {
          drawing = false;
          if((heightmap(j,i) > waterdepth || i == cols-1 ) && !na_matrix(j,i)) {
            if(i != cols-1) {
              double diff = heightmap(j,i)-heightmap(j,i-1);
              double adjustment_factor;
              if(diff == 0) {
                adjustment_factor = 0;
              } else {
                adjustment_factor = (waterdepth - heightmap(j,i-1))/diff;
              }
              endcoord = -(double)i - adjustment_factor;
            } else {
              endcoord = -cols;
            }
          } else {
            if(i != cols-1) {
              endcoord = -(double)i;
            } else {
              endcoord = -cols;
            }
          }
          vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,waterdepth,waterdepth,-startcoord-1,endcoord),2,3));
        }
      }
      if(!drawing && (j == 0 || j == rows - 1)) {
        if(((heightmap(j,i) < waterdepth) || 
           (na_matrix(j,i-offsetside) && heightmap(j,i+offsetside2) < waterdepth)) &&
           ((j == 0 && !na_matrix(1, i)) || (j == rows - 1 && !na_matrix(rows - 2, i)))) {
          if(!na_matrix(j,i-offsetside)) {
            if(i != 0) {
              double diff = heightmap(j,i)-heightmap(j,i-1);
              double adjustment_factor;
              if(diff == 0) {
                adjustment_factor = 0;
              } else {
                adjustment_factor = (waterdepth - heightmap(j,i-1))/diff;
              }
              startcoord = ((double)i-1) + adjustment_factor;
            } else {
              startcoord = 0;
            }
          } else {
            if(i != 0) {
              if(na_matrix(j,i)) {
                startcoord = (double)i+1;
              } else {
                startcoord = (double)i;
              }
            } else {
              startcoord = 0;
            }
          }
          drawing = true;
        }
      }
      //Interior
      if(drawing && j != 0 && j != rows - 1) {
        //Finish drawing if not NA AND
          //the back left or back right entries are NA AND the front left and front right entries are NOT NA OR
          //It is NA right in front of that entry
        if((!na_matrix(j,i) &&
           (((na_matrix(j+offset2,i-offsetside) || na_matrix(j-offset,i-offsetside)) && (!na_matrix(j+offset2,i+offsetside2) || !na_matrix(j-offset,i+offsetside2))) ||
           na_matrix(j,i+offsetside))) || i == cols - 1) {
          drawing = false;
          if(i != cols-1) {
            adjust = (waterdepth - heightmap(j,i-1))/(heightmap(j,i)-heightmap(j,i-1));
            if(heightmap(j,i) > waterdepth && fabs(adjust) < 1) {
              endcoord = (double)i + adjust;
            } else {
              endcoord = (double)i+1;
            }
          } else {
            endcoord = cols;
          }
          vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,waterdepth,waterdepth,-startcoord-1,-endcoord),2,3));
        }
      }
      if(!drawing && j != 0 && j != rows - 1) {
        //Start drawing if under water or the next space is underwater, but only if:
        //The matrix is not NA in the next entry AND
        //the current entry is not NA AND
        //the left OR right OR left front OR right front is NA
        if((heightmap(j,i) < waterdepth || (heightmap(j,i) >= waterdepth && heightmap(j,i+offsetside2) < waterdepth)) &&
           !na_matrix(j,i+offsetside2) &&
           (!na_matrix(j,i) && (na_matrix(j-offset,i) || na_matrix(j+offset2,i) || na_matrix(j+offset2,i+offsetside2) || na_matrix(j-offset,i+offsetside2)))) {
          if(i != 0) {
            adjust = (waterdepth - heightmap(j,i-1))/(heightmap(j,i)-heightmap(j,i-1));
            if(heightmap(j,i) > waterdepth && fabs(adjust) < 1) {
              startcoord = ((double)i-1) + adjust;
            } else {
              startcoord = (double)i;
            }
          } else {
            startcoord = 0;
          }
          drawing = true;
        }
      }
    }
  }
  for(int j = 0; j < cols; j++) {
    drawing = false;
    if(j != 0) {
      offset = 1;
    } else {
      offset = 0;
    }
    if(j != cols-1) {
      offset2 = 1;
    } else {
      offset2 = 0;
    }
    for(int i = 0; i < rows; i++) {
      if(i != 0) {
        offsetside = 1;
      } else {
        offsetside = 0;
      }
      if(i != rows-1) {
        offsetside2 = 1;
      } else {
        offsetside2 = 0;
      }
      //Edges
      if(drawing && (j == 0 || j == cols - 1)) {
        if(heightmap(i,j)  > waterdepth || i == rows-1 || na_matrix(i,j)) {
          drawing = false;
          if((heightmap(i,j)  > waterdepth || i == rows-1) && !na_matrix(i,j)) {
            if(i != rows-1) {
              double diff = heightmap(i,j)-heightmap(i-1,j);
              double adjustment_factor;
              if(diff == 0) {
                adjustment_factor = 0;
              } else {
                adjustment_factor = (waterdepth - heightmap(i-1,j))/diff;
              }
              endcoord = (double)i + adjustment_factor;
            } else {
              endcoord = rows;
            }
          } else {
            if(i != rows-1) {
              endcoord = (double)i;
            } else {
              endcoord = rows;
            }
          }
          vertices.push_back(vec2matrix(NumericVector::create(startcoord+1,endcoord,waterdepth,waterdepth,-1-j,-1-j),2,3));
        }
      }
      if(!drawing && (j == 0 || j == cols - 1)) {
        if((heightmap(i,j) < waterdepth || 
           (na_matrix(i - offsetside,j) && heightmap(i + offsetside2,j) < waterdepth)) &&
           ((j == 0 && !na_matrix(i, 1)) || (j == cols - 1 && !na_matrix(i, cols - 2)))) {
          if(!na_matrix(i-offsetside,j)) {
            if(i != 0) {
              double diff = heightmap(i,j)-heightmap(i-1,j);
              double adjustment_factor;
              if(diff == 0) {
                adjustment_factor = 0;
              } else {
                adjustment_factor = (waterdepth - heightmap(i-1,j))/diff;
              }
              startcoord = ((double)i-1) + adjustment_factor;
            } else {
              startcoord = 0;
            }
          } else {
            if(i != 0) {
              if(na_matrix(i,j)) {
                startcoord = (double)i+1;
              } else {
                startcoord = (double)i;
              }
            } else {
              startcoord = 0;
            }
          }
          drawing = true;
        }
      }
      //Interior
      if(drawing && j != 0 && j != cols - 1) {
        //Finish drawing if not NA AND
          //the back left or back right entries are NA AND the front left and front right entries are NOT NA OR
          //It is NA right in front of that entry
        if((!na_matrix(i,j) &&
           (((na_matrix(i-offsetside,j-offset) || na_matrix(i-offsetside,j+offset2)) && (!na_matrix(i+offsetside2,j-offset) || !na_matrix(i+offsetside2,j+offset2))) ||
           na_matrix(i+offsetside,j))) || i == rows - 1) {
          drawing = false;
          if(i != rows-1) {
            adjust = (waterdepth - heightmap(i-1,j))/(heightmap(i,j)-heightmap(i-1,j));
            if(heightmap(i,j) > waterdepth && fabs(adjust) < 1) {
              endcoord = (double)i + adjust;
            } else {
              endcoord = (double)i+1;
            }
          } else {
            endcoord = rows;
          }
          vertices.push_back(vec2matrix(NumericVector::create(startcoord+1,endcoord,waterdepth,waterdepth,-1-j,-1-j),2,3));
        }
      }
      if(!drawing && j != 0 && j != cols - 1) {
        //Start drawing if under water or the next space is underwater, but only if:
        //The matrix is not NA in the next entry AND
        //the current entry is not NA AND
        //the left OR right OR left front OR right front is NA
        if((heightmap(i,j) < waterdepth || (heightmap(i,j) >= waterdepth && heightmap(i+offsetside2,j) < waterdepth) ) && //Check depths
           !na_matrix(i+offsetside2,j) && //Not NA in the next entry
           (!na_matrix(i,j) && (na_matrix(i,j-offset) || na_matrix(i,j+offset2) || na_matrix(i+offsetside2,j-offset) || na_matrix(i+offsetside2,j+offset2)))) {
          if(i != 0) {
            adjust = (waterdepth - heightmap(i-1,j))/(heightmap(i,j)-heightmap(i-1,j));
            if(heightmap(i,j) > waterdepth && fabs(adjust) < 1) {
              startcoord = (double)i-1 + adjust;
            } else {
              startcoord = (double)i;
            }
          } else {
            startcoord = 0;
          }
          drawing = true;
          }
        }
      }
    }
  List vectorlist = wrap(vertices);
  return(vectorlist);
}


// [[Rcpp::export]]
List make_baselines_cpp(NumericMatrix& heightmap,
                         LogicalMatrix& na_matrix,
                         double waterdepth) {
  std::vector<NumericMatrix> vertices;
  int rows = heightmap.nrow();
  int cols = heightmap.ncol();
  int offset, offset2 = 0;
  int offsetside, offsetside2 = 0;
  bool drawing = false;
  double startcoord, endcoord = 1;

  for(int j = 0; j < rows; j++) {
    drawing = false;
    if(j != 0) {
      offset = 1;
    } else {
      offset = 0;
    }
    if(j != rows-1) {
      offset2 = 1;
    } else {
      offset2 = 0;
    }
    for(int i = 0; i < cols; i++) {
      if(i != 0) {
        offsetside = 1;
      } else {
        offsetside = 0;
      }
      if(i != cols-1) {
        offsetside2 = 1;
      } else {
        offsetside2 = 0;
      }
      //Edges
      if(drawing && (j == 0 || j == rows - 1)) {
        if(i == cols-1 || na_matrix(j,i)) {
          drawing = false;
          if(i != cols-1) {
            endcoord = -(double)i;
          } else {
            endcoord = -cols;
          }
          vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,waterdepth,waterdepth,-startcoord-1,endcoord),2,3));
        }
      }
      if(!drawing && (j == 0 || j == rows - 1)) {
        if((i == 0 || na_matrix(j,i-offsetside)) && !na_matrix(j,i)) {
          if(i != 0) {
            startcoord = (double)i;
          } else {
            startcoord = 0;
          }
          drawing = true;
        }
      }
      //Interior
      if(drawing && j != 0 && j != rows - 1) {
        //Finish drawing if not NA AND
        //the back left or back right entries are NA AND the front left and front right entries are NOT NA OR
        //It is NA right in front of that entry
        if((!na_matrix(j,i) &&
           (((na_matrix(j+offset2,i-offsetside) || na_matrix(j-offset,i-offsetside)) && (!na_matrix(j+offset2,i+offsetside2) || !na_matrix(j-offset,i+offsetside2))) ||
           na_matrix(j,i+offsetside2))) || i == cols - 1) {
          drawing = false;
          if(i != cols-1) {
            endcoord = (double)i+1;
          } else {
            endcoord = cols;
          }
          vertices.push_back(vec2matrix(NumericVector::create(1+j,1+j,waterdepth,waterdepth,-startcoord-1,-endcoord),2,3));
        }
      }
      if(!drawing && j != 0 && j != rows - 1) {
        //Start drawing if under water or the next space is underwater, but only if:
        //The matrix is not NA in the next entry AND
        //the current entry is not NA AND
        //the left OR right OR left front OR right front is NA
        if((!na_matrix(j,i+offsetside2) &&
           (!na_matrix(j,i) && (na_matrix(j-offset,i) ||
           na_matrix(j+offset2,i) ||
           na_matrix(j+offset2,i+offsetside2) ||
           na_matrix(j-offset,i+offsetside2))))) {
          if(i != 0) {
            startcoord = (double)i;
          } else {
            startcoord = 0;
          }
          drawing = true;
        }
      }
    }
  }
  for(int j = 0; j < cols; j++) {
    drawing = false;
    if(j != 0) {
      offset = 1;
    } else {
      offset = 0;
    }
    if(j != cols-1) {
      offset2 = 1;
    } else {
      offset2 = 0;
    }
    for(int i = 0; i < rows; i++) {
      if(i != 0) {
        offsetside = 1;
      } else {
        offsetside = 0;
      }
      if(i != rows-1) {
        offsetside2 = 1;
      } else {
        offsetside2 = 0;
      }
      //Edges
      if(drawing && (j == 0 || j == cols - 1)) {
        if(i == rows-1 || na_matrix(i,j)) {
          drawing = false;
          if(i != rows-1) {
            endcoord = (double)i;
          } else {
            endcoord = rows;
          }
          vertices.push_back(vec2matrix(NumericVector::create(startcoord+1,endcoord,waterdepth,waterdepth,-1-j,-1-j),2,3));
        }
      }
      if(!drawing && (j == 0 || j == cols - 1)) {
        if((i == 0 || na_matrix(i-offsetside,j)) && !na_matrix(i,j)) {
          if(i != 0) {
            startcoord = (double)i;
          } else {
            startcoord = 0;
          }
          drawing = true;
        }
      }
      //Interior
      if(drawing && j != 0 && j != cols - 1) {
        //Finish drawing if not NA AND
        //the back left or back right entries are NA AND the front left and front right entries are NOT NA OR
        //It is NA right in front of that entry
        if((!na_matrix(i,j) &&
           (((na_matrix(i-offsetside,j-offset) || na_matrix(i-offsetside,j+offset2)) && (!na_matrix(i+offsetside2,j-offset) || !na_matrix(i+offsetside2,j+offset2))) ||
           na_matrix(i+offsetside2,j)))  || i == rows - 1) {
          drawing = false;
          if(i != rows-1) {
            endcoord = (double)i+1;
          } else {
            endcoord = rows;
          }
          vertices.push_back(vec2matrix(NumericVector::create(startcoord+1,endcoord,waterdepth,waterdepth,-1-j,-1-j),2,3));
        }
      }
      if(!drawing && j != 0 && j != cols - 1) {
        //Start drawing if under water or the next space is underwater, but only if:
        //The matrix is not NA in the next entry AND
        //the current entry is not NA AND
        //the left OR right OR left front OR right front is NA
        if(!na_matrix(i+offsetside2,j) &&
           (!na_matrix(i,j) && (na_matrix(i,j-offset) || na_matrix(i,j+offset2) || na_matrix(i+offsetside2,j-offset) || na_matrix(i+offsetside2,j+offset2)))) {
          if(i != 0) {
            startcoord = (double)i;
          } else {
            startcoord = 0;
          }
          drawing = true;
        }
      }
    }
  }
  List vectorlist = wrap(vertices);
  return(vectorlist);
}
