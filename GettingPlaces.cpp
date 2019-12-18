/****************************************************************************************
                        	       Getting Places
                            	     by Luis Sosa

****************************************************************************************/

#include "library.h"
#include <limits>
static const int INF = (numeric_limits<int>::max)();

//////////////////////////////////////////
/// HASHTABLE ///////////////////////////
/////////////////////////////////////////

int hashify(string s){
  //Hashing constants
  const long long H0 = 98297;
  const long long H1 = 245209;
  int h = H0;

  for(int i = 0; i<s.length(); ++i){
    h = h*H1 + s[i];
  }
  if(h<0){
    h = -h;
  }
 // cout << h << "\n";
  return h;
}

void trimSpaces (string &s){
  int length = s.length();

  //Trim right end spaces
  for (int i=length-1; i>=0; --i){
    if (s[i] == ' '){
      s.erase(i,1);
      --length;
    }
    else{
     break;
    }
  }

  //Trim left end spaces
  while(length>0){
    if (s[0] == ' '){
      s.erase(0,1);
      --length;
    }
    else{
      break;
    }
  }
}

struct Place {
  int id;
  string state;
  string name;
  int population;
  double area;
  double latitude;
  double longitude;
  int intersectionID;
  double distanceToRoad;
  Place (int i, string s, string n, int p, double a, double la, double lo, int in, double d){
	id = i;
    state = s;
    name = n;
    population = p;
    area = a;
    latitude = la;
    longitude = lo;
    intersectionID = in;
    distanceToRoad = d;
  }

  void print(){
    cout << id << " " << state << " " << name << " " << population << " "
         << area << " " << latitude << " " << longitude << " " << intersectionID
         << " " << distanceToRoad << "\n";
  }

};

struct HashLink{
  string key;
  HashLink* next;
  Place* place;
  HashLink(Place* p, HashLink* n = NULL){
    next = n;
    place = p;
    key = p->name + p->state;
  }
};

struct HashTable{
  HashLink** table;
  int tableSize;

  HashTable(const int s){
    tableSize = s;
    table = new HashLink*[s];
    for(int i=0; i<s; ++i){
      table[i] = NULL;
    }
  }

  ~HashTable(){
    for(int i=0; i<tableSize; ++i){
      HashLink* current = table[i];
      while(current != NULL){
        HashLink* following = current->next;
        delete current;
        current = following;
      }
    }
    delete[] table;
  }

  void addToTable(HashLink* hl){
    int position = hashify(hl->key) % tableSize;
    HashLink* firstLink = table[position];
    hl->next = firstLink;
    table[position] = hl;
  }

  HashLink* find(string name, string state){
    string query = name + state;
    int location = hashify(query) % tableSize;
    HashLink* current = table[location];
    while(current != NULL){
      if(current->key == query){
        return current;
      }
      else{
        current = current->next;
      }
    }
    return NULL;
  }

  HashLink* operator[](int pos){
    return table[pos];
  }
};

int getFileSize(ifstream &fin){
  fin.seekg(0,fin.end);
  int size = fin.tellg();
  fin.seekg(0, fin.beg);
  return size;
}

string getParameter(string &line, const int s[], int &i){
  string parameter = line.substr (s[i],s[i+1]);
  trimSpaces(parameter);
  i+=2;
  return parameter;
}

void readNamedPlaces(HashTable &hashTable,const char* fileName,const int s[], const int HASHSIZE){
  ifstream fin;
  fin.open(fileName);
  if(!fin){
    cerr << "Error opening file\n";
    exit(-1);
  }
  int fileSize = getFileSize(fin);

  string line;
  int id = 0;
  while(getline(fin, line)){
	  int i = 0;

    string tempID = getParameter(line, s, i);
    int id = atoi(tempID.c_str());

    string state = getParameter(line, s, i);
	string name = getParameter(line, s, i);
    string tempPop = getParameter(line, s, i);
    int pop = atoi(tempPop.c_str());

    string tempArea = getParameter(line, s, i);
    double area = atof(tempArea.c_str());

    string tempLati = getParameter(line, s, i);
    double lati = atof(tempLati.c_str());

    string tempLongi = getParameter(line, s, i);
    double longi = atof(tempLongi.c_str());

    string tempInter = getParameter(line, s, i);
    int inter = atoi(tempInter.c_str());

    string tempDist = getParameter(line, s, i);
    double dist = atof(tempDist.c_str());

    Place* p = new Place(id,state,name,pop,area,lati,longi,inter,dist);
    HashLink* hl = new HashLink(p);
    hashTable.addToTable(hl);
  }

  fin.close();
}

//////////////////////////////////////////////////
// SHORTEST PATH ////////////////////////////////
////////////////////////////////////////////////

struct Road;
struct Intersection;

struct Road{
  string name;
  int length;
  int endA;
  int endB;

  Road(string n, double l, int a, int b){
    name = n;
    length = l*1000; //Internally stored as an int for easy comparisons
    endA = a;
    endB = b;
  }

};

struct Intersection{
  int id;
  string name;
  string state;
  double distance;
  double latitude;
  double longitude;
  int numRoads;
  vector<Road*> roads;
  int priority;
  int heapPos;

  Intersection(int i, string n, string s, double lo, double la, double d){
    id = i;
    name = n;
    state = s;
    distance = d;
    longitude = lo;
    latitude = la;
    numRoads = 0;
    priority = INF;
    heapPos = -1; //-1: not explored, -2: already explored
  }

  void addRoad(Road* r){
      roads.push_back(r);
      ++numRoads;
  }

  void setPriority(int p){
    priority = p;
  }

  void setHeapPos(int p){
    heapPos = p;
  }

  void reset(){
	  priority = INF;
	  heapPos = -1;
  }

};

void swapI(Intersection* &a, Intersection* &b){
  int oldAHeapPos = a->heapPos;
  a->heapPos = b->heapPos;
  b->heapPos = oldAHeapPos;
  
  Intersection* oldA = a;
  a = b;
  b = oldA;
}

struct Map{
  vector<Intersection*> intersections;
  vector<Road*> roads;
  
  void add(Intersection* i){
    intersections.push_back(i);
  }

  void add(Road* r){
    roads.push_back(r);
    intersections[r->endA]->roads.push_back(r);
    ++intersections[r->endA]->numRoads;
    intersections[r->endB]->roads.push_back(r);
    ++intersections[r->endB]->numRoads;
  }

  void reset(){
	  for(int i=0; i<intersections.size(); ++i){
		  intersections[i]->reset();
	  }
  }
};



typedef struct MinHeap{
  vector<Intersection*> intersections;
  int size;

  MinHeap(){
    size = 0;
  }

  // Adjusts position of inserted/updated entry
  // returns final position
  int insertAdjust(int current){
    // If entry is at the top of the heap, that is its final position
    if(current <= 0){
      return current;
    }

    // Make entries float if priority is lower than parents
    int parent = (current-1)/2;
    if(intersections[current]->priority < intersections[parent]->priority){
      swapI(intersections[current], intersections[parent]);
      return insertAdjust(parent);
    }

    // Else current position is final
    return current;
  }

  // Priority should only decrease
  // Returns true if priority was changed, false if it was not
  bool changePriority(int heapPos, int newPriority){
    if(newPriority < intersections[heapPos]->priority){
      intersections[heapPos]->priority = newPriority;
      return true;
    }
    else return false;
  }

  void add(Intersection* i, int priority){
    // If place is already in the heap
    if(i->heapPos >= 0){
      if(changePriority(i->heapPos, priority)){
        insertAdjust(i->heapPos);
      }
    }
    // If place is not yet in the heap
    else if(i->heapPos == -1){
      i->setPriority(priority);
      intersections.push_back(i);
      i->setHeapPos(size);
      insertAdjust(size);
      ++size;
    }
  }

  // Adjusts heap after removal of first entry
  void popAdjust(int current = 0){
    int leftChild = current*2 + 1;
    if(leftChild >= size){
      return;
    }
    int rightChild = current*2 + 2;
    if(rightChild >= size){
      rightChild = -1;
    }

    // If there is no right child
    if(rightChild < 0){
      if(intersections[leftChild]->priority < intersections[current]->priority){
        swapI(intersections[leftChild],intersections[current]);
        return;
      }
    }

    else{
      // Swap with left child
      if(intersections[leftChild]->priority < intersections[current]->priority
        && intersections[leftChild]->priority <= intersections[rightChild]->priority){
          swapI(intersections[leftChild], intersections[current]);
          popAdjust(leftChild);
      }
      // Swap with right child
      else if(intersections[rightChild]->priority < intersections[current]->priority
        && intersections[rightChild]->priority <= intersections[leftChild]->priority){
           swapI(intersections[rightChild], intersections[current]);
           popAdjust(rightChild);
      }
    }
  }

  // Remove element with least priority
  Intersection* pop(){
    if(size<=0){
      return NULL;
    }
    --size;
    swapI(intersections[size], intersections[0]);
    popAdjust();
    Intersection* popped = intersections.back();
    popped->heapPos = -2;
    intersections.pop_back();
    return popped;
  }

  Intersection* operator[](int i){
    return intersections[i];
  }

} priorityQueue;



// Create intersection objects from file and add them to map
void readIntersections(Map &map, const char* filename){
  ifstream fin;
  fin.open(filename);
  if(fin.fail()){
    cerr << "Error opening intersections file";
    exit(-1);
  }
  string line;
  int id = 0;
  while(getline(fin, line)){

    string d = line.substr(19,8);
    trimSpaces(d);
    double dist = atof(d.c_str());

    string lo = line.substr(0,10);
    trimSpaces(lo);
    double longi = atof(lo.c_str());

    string la = line.substr(10,10);
    trimSpaces(la);
    double lati = atof(la.c_str());

    string state = line.substr(28,2);

    string name = line.substr(30,line.length()-1);
    trimSpaces(name);

    Intersection* i = new Intersection(id, name, state, longi, lati, dist);
    //cout << "Location " << id << ": " << i->distance << " miles away from " << i->name << ", " << i->state << "\n";
    map.add(i);
    ++id;
  }
  fin.close();
}

// Create connection objects from file and add them to map
void readConnections(Map &map, const char* filename){
  ifstream fin;
  fin.open(filename);
  if(fin.fail()){
    cerr << "Error opening connections file";
    exit(-1);
  }
  string name;
  while(fin >> name){
    string code;
    int endA;
    int endB;
    double length;
    fin >> code >> endA >> endB >> length;

    if(name == "?"){
      name = "unnamed";
      if(code[0] != 'F'){
        name += " road";
      }
    }
    else if(code[0] != 'F'){
      name = "road " + name;
    }

    if(code[0] == 'F'){
      name += " ferry crossing";
    }

    Road* r = new Road(name, length, endA, endB);
    //cout << r->name << " of length " << r->length << " connecting " << r->endA << " and " << r->endB << "\n";
    map.add(r);
  }
  fin.close();
}

// Given a location and a road, output the other location
int otherEnd(Intersection* current, Road* road){
  if(road->endA == current->id){
    return road->endB;
  }
  else{
    return road->endA;
  }
}

// Continues exploration, returns shortest distance
int findShortestPath(Map &map, priorityQueue &pq, Intersection* end){
  Intersection* current = pq.pop();

  //If no route was found
  if(current == NULL){
    return -1;
  }

  //cerr << "Current location: " << current->id << " ";
  int totalDist = current->priority;
  //cerr << " d = " << totalDist << "\n";

  //If destination was reached
  if(current->id == end->id){
    return totalDist;
  }

  //Add all contiguous locations to the priority queue
  for(int i=0; i < current->numRoads; ++i){
    int nextID = otherEnd(current, current->roads[i]);
    Intersection* next = map.intersections[nextID];
    pq.add(next, totalDist + current->roads[i]->length);
  }

  return findShortestPath(map, pq, end);
}

// Initiates exploration
int findShortestPath(Map &map, priorityQueue &pq, Intersection* start, Intersection* end){
  pq.add(start, 0);
  return findShortestPath(map, pq, end);
}

// Checks if a segment belongs to the shortest path, taking into account rounding errors
bool isShortestPath(int i, int distance, Intersection* current, Intersection* next){
  //cout << "This Priority " << distance << ", NextPriority " << next->priority << ", road dist: " << current->roads[i]->length << "\n";
  return (distance - current->roads[i]->length == next->priority);
}

// Returns approximate cardinal direction between two points
string getDirection(Intersection* from, Intersection* to){
  double x1 = from->longitude;
  double y1 = from->latitude;
  double x2 = to->longitude;
  double y2 = to->latitude;

  //Determine direction based on slope of the line that connects "from" and "to"
  double m = (y2-y1) / (x2-x1);
  
  // Right semicircle
  if(x2 > x1){
    if(m > 2) return "N";
    else if(m<=2 && m>=0.5) return "NE";
    else if(m<0.5 && m>-0.5) return "E";
    else if(m>=-2 && m<=-0.5) return "SE";
    else return "S";
  }
  
  // Left cemicircle
  else if(x2 < x1){
    if(m > 2) return "S";
    else if(m<=2 && m>=0.5) return "SW";
    else if(m<0.5 && m>-0.5) return "W";
    else if(m>=-2 && m<=-0.5) return "NW";
    else return "N";
  }
  
  //If the points lie on the same vertical, distinguish by y values
  else if(y2 > y1){
    return "N";
  }
  else{
    return "S";
  }
}

////////////////////////
///// Graphics /////////
///////////////////////

string scanCoverage(Place* endpoints[], string filename="./MapInfo/coverage.txt"){
	double lati1 = endpoints[0]->latitude;
	double lati2 = endpoints[1]->latitude;
	double longi1 = endpoints[0]->longitude;
	double longi2 = endpoints[1]->longitude;

	ifstream fin(filename);
	if(fin.fail()){
		cerr << "Error opening coverage file " << filename << "\n"; 
		exit(-1);
	}
	int maxLati, minLati, minLongi, maxLongi;
	string file;
	string smallestFile = "";
	int smallestSize = INF;
	while(fin >> maxLati){
		fin >> minLati >> minLongi >> maxLongi >> file;
		if(fin.fail()){
			cerr << "Error reading file " << filename << "\n";
			exit(-2);
		}
		
		if(lati1 >= minLati && lati1 <= maxLati && lati2 >= minLati && lati2 <= maxLati &&
			longi1 >= minLongi && longi1 <= maxLongi && longi2 >= minLongi && longi2 <= maxLongi){
			int size = (maxLati-minLati)*(maxLongi-minLongi);
			if( size < smallestSize){
				smallestFile = file;
				smallestSize = size;
			}
		}
	}
	fin.close();
	return "./Tiles/" + smallestFile;
}

struct MapInfo{
	int numRows;
	int numCols;
	int pixelBytes;
	int pixelSeconds;
	int leftmostSecond;
	int topSecond;
	short minAltitude;
	short maxAltitude;
	short oceanValue;

	MapInfo(int r, int c, int pb, int ps, int l, int t, short m, short M, short o){
		numRows = r;
		numCols = c;
		pixelBytes = pb;
		pixelSeconds = ps;
		leftmostSecond = l;
		topSecond = t;
		minAltitude = m;
		maxAltitude = M;
		oceanValue = o;
	}
};

void drawMapPixel(short altitude, MapInfo* info){
	if(altitude == info->oceanValue){
		if(random_in_range(0,10))
			set_pixel_color(color::dark_blue);
		else
			set_pixel_color(color::light_blue);
	}
	else{
		if(altitude < 2)
			set_pixel_color(make_color_int(107, 255, 179));
		else if(altitude < 50)
			set_pixel_color(make_color_int(0,230,0));
		else if(altitude < 75)
			set_pixel_color(make_color_int(153,255,102));
		else if(altitude < 140)
			set_pixel_color(make_color_int(0, 204, 0));
		else if(altitude < 300)
			set_pixel_color(make_color_int(76, 230, 0));
		else if(altitude < 500)
			set_pixel_color(make_color_int(102, 204, 0));
		else if(altitude < 700)
			set_pixel_color(make_color_int(255, 219, 77));
		else if(altitude < 1000)
			set_pixel_color(make_color_int(191, 128, 64));
		else if(altitude < 1500)
			set_pixel_color(make_color_int(228,205,108));
		else if(altitude < 1700)
			set_pixel_color(make_color_int(137,91,42));
		else if(altitude < 1900)
			set_pixel_color(color::light_grey);
		else
			set_pixel_color(color::white);
	}
}

MapInfo* drawMap(string tileFilename){
	ifstream fin(tileFilename, ios::binary);

	if(fin.fail()){
		cerr << "Error opening maptile file\n";
		exit(-4);
	}

	//Get file metadata
	int numRows;
	int numCols;
	int pixelBytes;
	int pixelSeconds;
	int leftmostSecond;
	int topSecond;
	short minAltitude;
	short maxAltitude;
	short oceanValue;

	string ignore;
	fin >> ignore >> numRows
	    >> ignore >> numCols
		>> ignore >> pixelBytes
		>> ignore >> pixelSeconds
		>> ignore >> leftmostSecond
		>> ignore >> topSecond
		>> ignore >> minAltitude
		>> ignore >> maxAltitude
		>> ignore >> oceanValue;

	MapInfo* info = new MapInfo(numRows, numCols, pixelBytes, pixelSeconds, leftmostSecond, topSecond, minAltitude, maxAltitude, oceanValue);

	// Draw map as it is being read
	make_window(numCols-1, numRows-1);
	short altitude;
	fin.seekg(numCols * sizeof(short));
	for(int r=0; r<numRows; ++r){
		for(int c=0; c<numCols; ++c){
			fin.read((char*)&altitude,sizeof(short));
			move_to(c,r);
			drawMapPixel(altitude, info);
			
		}
	}
	fin.close();
	return info;
}

int longitudeToPixels(double longi, MapInfo* info){
	return (3600*longi - info->leftmostSecond)/info->pixelSeconds;
}

int latitudeToPixels(double lati, MapInfo* info){
	return (info->topSecond - 3600*lati)/info->pixelSeconds;
}

void drawRoute(Intersection* current, Road* previous, MapInfo* info, bool final=false){
	int y = latitudeToPixels(current->latitude, info);
	int x = longitudeToPixels(current->longitude, info);

	if(previous == NULL){
		draw_to(x,y);
		set_pen_width(12);
		set_pen_color(color::black);
		move_to(x,y);
		draw_point();
		set_pen_width(10);
		set_pen_color(color::red);
		draw_point();
	}
	else if(final){
		set_pen_width(12);
		set_pen_color(color::black);
		draw_point(x,y);
		set_pen_width(10);
		set_pen_color(color::light_blue);
		draw_point(x,y);
		set_pen_width(2);
		set_pen_color(color::black);
		draw_to(x,y);
	}
	else{
		draw_to(x,y);
	}
}

// Prints navigation instructions
// Since the shortest path is inferred in the inverse order it is travelled, previous holds the next position to travel to
void navigate(Map &map, MapInfo* info, Intersection* start, Intersection* current, int distance, Intersection* previous = NULL, Road* segment = NULL, double segmentLength = 0){

  // Base case -- Starting location reached
  if(distance <= 0){
	drawRoute(current, segment, info, true);
    cout << "Start by heading " << segmentLength/1000.0 << " miles " << getDirection(current,previous)
         << " from " << current->name << " to " << previous->name << " via " << segment->name << "\n";
    return;
  }

  // Obtain segments of shortest path
  for(int i=0; i<current->numRoads; ++i){
    Intersection* next = map.intersections[otherEnd(current, current->roads[i])];

   // cerr << distance-current->roads[i]->length << " looking for\n";
   // cerr << next->priority << " looking at\n";

    if(isShortestPath(i, distance, current, next)){
      distance-= current->roads[i]->length;

	  // New Segment
	  if(segment != NULL && current->roads[i]->name != segment->name){
        navigate(map, info, start, next, distance, current, current->roads[i], current->roads[i]->length);
		
		if(previous != NULL){
		  cout << "Next, continue by heading " << segmentLength/1000.0 << " miles " << getDirection(current,previous)
               << " from " << current->name << " to " << previous->name << " via " << segment->name << "\n";
		}
	  }

	  // Continue Previous Segment
	  else{
		segmentLength += current->roads[i]->length;
		navigate(map, info, start, next, distance, previous, current->roads[i], segmentLength);
	  }
     break;
    }
  }

  // Final print message
  if(segment == NULL){
    cout << "Congratulations! You have arrived at your destination: " << current->name << "\n";
  }

  drawRoute(current, segment, info, false);
}

void setTitle(string start, string end){
		string title = "From " + start + " to " + end;
		set_caption(title);
}


void main(){

	cout << "Booting up. Please wait a few seconds...\n";

	// Load Named Places
	const int HASHSIZE = 4096;
	HashTable hashTable(HASHSIZE);
	const int COLUMNSIZES[] = {0,8, 8,2, 10,49, 59,7, 66,14, 80,9, 90,11, 101,5, 106,8};

	string  fileName = "./MapInfo/named-places.txt";
	readNamedPlaces(hashTable, fileName.c_str(), COLUMNSIZES, HASHSIZE);

	// Load Intersections and Connections
	Map map;

	string interFileName = "./MapInfo/intersections.txt";
	readIntersections(map, interFileName.c_str());

	string connecFileName = "./MapInfo/connections.txt";
	readConnections(map, connecFileName.c_str());

	
	
	
begin:
	while(true){
		priorityQueue pq;

		// User Input
		string trip[4]; //startName, startState, endName, endState
		cout << "\nPlease input a starting name and state, and a destination name and state\n" 
		<< "Please separate your four inputs using commas. Eg. Miami Beach, FL, New York, NY\n"
		<< "Click on the upper right X to quit\n";

		for(int i=0; i<3; ++i){
			 getline(cin, trip[i], ',');
			 trimSpaces(trip[i]);
		}
		cin >> trip[3];
		cin.ignore(10000,'\n');

		// Determine trip endpoints
		Place* endpoints[2];
		int i = 0;
		for(int j=0; j<4; j+=2){
			HashLink* result = hashTable.find(trip[j],trip[j+1]);
			if(result == NULL){
				cerr << "The place " << trip[j] << ", " << trip[j+1] << " could not be found\n";
				goto begin;
			}
			else{
				endpoints[i] = result->place;
				++i;
			}
		}
		// Determine map file to be used
		string tileFile = scanCoverage(endpoints);
		if(tileFile == ""){
			cerr << "Could not find a map that contains both " <<  trip[0] << ", " << trip[1] << " and "
				 << trip[2] << ", " << trip[3] << "\n";
			continue;
		}

		// Draw selected map file
		MapInfo* info = drawMap(tileFile);

		// Find shortest path
		Intersection* start =  map.intersections[endpoints[0]->intersectionID];
		Intersection* end = map.intersections[endpoints[1]->intersectionID];
		setTitle(start->name, end->name);

		int dist = findShortestPath(map, pq, start, end);
		if(dist < 0){
			cout << "There is no known way to get from " << trip[0] << ", " << trip[1] << " to "
				 << trip[2] << ", " << trip[3] << "\n";
			continue;
		}

		// Draw Route
		cout << "Your total trip distance will be " << dist/1000.0 << " miles ";
			if(dist>20)
				cout << "or roughly " << dist/64000 << " hours, " << floor((dist/64000.0 - dist/64000)*60) << " minutes\n";
			else
				cout << "\n";
		navigate(map, info, start, end, dist);
		map.reset();
	}
}


