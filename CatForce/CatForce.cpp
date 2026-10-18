// CatForce - Catalyst search utility based on LifeAPI using brute force.
// Written by Michael Simkin 2015
#include "LifeAPI.h"
#include "RuleCompiler.h"
#include <unistd.h>
#include <sys/stat.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <csignal>
#include <unordered_map>
#include <map>
#include <set>
#include <array>
#include <algorithm>

const int MAX_SEARCHED = 5;   // num-catalyst limit
const int MAX_CATALYSTS = 8;  // searched + fixed-cat
const int MAX_BAITS = 8;
const unsigned kFilterMaxGen = 1999999999u; // filter range ending in "max", until max-gen is known

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
}

std::vector<std::string> splitwhitespace(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
                                  std::istream_iterator<std::string>());
    return ret;
}

enum FilterType {
  ANDFILTER,
  ORFILTER,
  MATCHFILTER,
};

class SearchParams {
public:
  unsigned maxGen;
  unsigned numCatalysts;      // the most catalysts per result (num-catalyst N or M-N)
  unsigned minCatalysts = 0;  // the fewest (0 = same as numCatalysts)
  std::string sortBy = "found"; // output order, see SortResults
  unsigned numTransparent;
  unsigned stableInterval;
  std::string pat;
  int xPat;
  int yPat;
  unsigned startGen;
  unsigned lastGen;
  std::string outputFile;
  std::string fullReportFile;
  int searchArea[4]{};
  // Search area = union of the rectangles in searchAreas (or the default
  // searchArea if none were given) minus the union of excludeAreas.
  std::vector<std::array<int, 4>> searchAreas;
  std::vector<std::array<int, 4>> excludeAreas;
  int maxW;
  int maxH;
  StaticSymmetry symmetry;
  std::vector<SymmetryTransform> symmetryChain;
  std::string symmetryName = "C1";
  std::vector<std::string> targetFilter;
  std::vector<int> filterdx;
  std::vector<int> filterdy;
  std::vector<int> filterGen;
  std::vector<std::pair<int, int>> filterGenRange;
  std::vector<FilterType> filterType;

  int maxCatSize;

  std::string alsoRequired;
  std::pair<int, int> alsoRequiredXY;

  // "bait": an object in the starting state that must be destroyed/disturbed
  // at some point and later reappear and survive.
  struct BaitInput {
    std::string rle;
    int x = 0, y = 0;
    int survive = -1;  // -1: use stable-interval
    int absent = 1;    // minimum number of generations it must be missing
  };
  std::vector<BaitInput> baits;

  int stopAfterCatsDestroyed;
  int maxJunk;
  bool keepDuplicates;
  bool keepInteracting = false; // keep solutions whose catalysts interact with each other
  int outputGap;
  unsigned numThreads = 0;   // 0 = one per hardware thread (default)
  bool verbose = false;      // print every top-level placement (the old behaviour)
  unsigned partIndex = 0, numParts = 1; // run only part K of N of the search
  std::string inputFile;
  unsigned maxTimeSec = 0;   // 0 = unlimited
  unsigned maxResults = 0;   // 0 = unlimited
  std::string outputList;    // optional text listing of every solution
  unsigned long resumeFrom = 0; // skip this many top-level placements (from a paused run)
  unsigned splitLevel = 1;   // 2 = threads share out second-level placements
  bool soloRecovery = false; // auto contact filter: a new catalyst must recover alone at least once
  bool lightcone = true;
  bool lightconeTight = false; // unsafe: failure region = only the mismatched cells at the failure time     // prune placements that cannot reach the next catalyst failure in time
  std::string uniqueFile;    // one square per distinct outcome (default: <output>_unique.rle)
  // Catalysts are not placed where their first contact would be with a spaceship
  // that is more than this many cells away from every other non-catalyst cell
  // (an escaped spaceship, e.g. a glider being eaten far from the reaction). 0 = off.
  int spaceshipCutoff = 10;
  int torus = N;             // torus size; a different value re-runs the CatForce<size> binary
  // An RLE placed with its top-left corner at (x, y)
  struct PlacedRle {
    std::string rle;
    int x = 0, y = 0;
    std::string rest; // fixed-cat: max-active and catalyst attributes
  };
  std::vector<PlacedRle> fixedCats;      // catalysts already in place (not searched for)
  std::vector<PlacedRle> catalystZones;  // catalysts must lie inside these cells
  std::vector<PlacedRle> noCatalystZones; // ... and outside these
  std::vector<PlacedRle> keepEmpty;      // cells that must stay dead in every generation
  std::string torusEdge = "flag";        // results reaching the torus edge: flag | drop | off
  bool patSeen = false;

  SearchParams() {
    maxGen = 250;
    numCatalysts = 2;
    numTransparent = 100;
    stableInterval = 15;
    pat = "";
    searchArea[0] = -30;
    searchArea[1] = -30;
    searchArea[2] = 60;
    searchArea[3] = 60;
    xPat = 0;
    yPat = 0;
    startGen = 1;
    lastGen = 100;
    outputFile = "results.rle";
    fullReportFile = "";
    maxW = -1;
    maxH = -1;
    symmetry = StaticSymmetry::C1;
    symmetryChain = {};
    maxCatSize = -1;
    alsoRequired = "";
    alsoRequiredXY = {0, 0};
    stopAfterCatsDestroyed = -1;
    maxJunk = -1;
    keepDuplicates = false;
    outputGap = 10;
  }
};

// ---- input validation helpers -------------------------------------------
[[noreturn]] static void Die(const std::string &msg) {
  std::cout << msg << std::endl;
  exit(1);
}

// Whole number with a minimum value; atoi() silently turned typos into 0.
static int ParseIntOrDie(const std::string &s, const std::string &what, long minVal = -1000000000L) {
  char *end = nullptr;
  long v = strtol(s.c_str(), &end, 10);
  if (s.empty() || *end != '\0')
    Die(what + ": \"" + s + "\" is not a whole number");
  if (v < minVal)
    Die(what + ": " + s + " is below the minimum of " + std::to_string(minVal));
  if (v > 1000000000L)
    Die(what + ": " + s + " is too large");
  return (int)v;
}

// "90", "90s", "15m", "2h" or "1d" -> seconds
static unsigned ParseDuration(const std::string &str, const std::string &what) {
  std::string num = str;
  unsigned mult = 1;
  if (!num.empty() && isalpha((unsigned char)num.back())) {
    char u = (char)tolower((unsigned char)num.back());
    if (u == 's') mult = 1;
    else if (u == 'm') mult = 60;
    else if (u == 'h') mult = 3600;
    else if (u == 'd') mult = 86400;
    else Die(what + ": unknown unit in \"" + str + "\" (use s, m, h or d)");
    num.pop_back();
  }
  return (unsigned)ParseIntOrDie(num, what, 1) * mult;
}

// "2/4" -> run the 2nd of 4 parts (stored 0-based)
static void ParsePart(const std::string &str, unsigned &index, unsigned &count) {
  size_t slash = str.find('/');
  if (slash == std::string::npos)
    Die("part: \"" + str + "\" must look like K/N, e.g. 1/4");
  int k = ParseIntOrDie(str.substr(0, slash), "part K", 1);
  int n = ParseIntOrDie(str.substr(slash + 1), "part N", 1);
  if (k > n)
    Die("part: K must be between 1 and N (got " + str + ")");
  index = (unsigned)(k - 1);
  count = (unsigned)n;
}

// Inserts a suffix before a file name's extension: a/b.rle + _x -> a/b_x.rle
static std::string AddSuffix(const std::string &name, const std::string &suffix) {
  size_t slash = name.find_last_of("/\\");
  size_t dot = name.find_last_of('.');
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash))
    return name + suffix;
  return name.substr(0, dot) + suffix + name.substr(dot);
}

// Writes a text file via a temporary name + rename, so a crash or a second
// Ctrl-C during a save never leaves a half-written results file.
static void SaveTextAtomically(const std::string &fname, const std::string &text) {
  std::string tmp = fname + ".saving";
  {
    std::ofstream out(tmp.c_str(), std::ofstream::binary);
    out << text;
    if (!out.good()) {
      out.close();
      std::remove(tmp.c_str());
      std::cout << "Could not write " << fname << std::endl;
      return;
    }
  }
  if (std::rename(tmp.c_str(), fname.c_str()) != 0) {
    std::remove(tmp.c_str());
    std::cout << "Could not write " << fname << std::endl;
  }
}

// The RLE parser only understands digits, o, b, $ and !. Anything else used to
// be silently dropped, which also skipped the pattern's x/y offset.
static void CheckRle(const std::string &rle, const std::string &what) {
  if (rle.empty())
    Die(what + ": empty RLE");
  for (char c : rle)
    if (!((c >= '0' && c <= '9') || c == 'o' || c == 'b' || c == '$' || c == '!'))
      Die(what + ": unsupported character '" + std::string(1, c) + "' in RLE \"" + rle +
          "\" (only o, b, $, digits and ! are allowed, with no spaces or x = ... header)");
}

class CatalystInput {
public:
  std::string rle;
  unsigned maxDisappear;
  int centerX;
  int centerY;
  char symmType;
  std::vector<std::string> forbiddenRLE;
  std::vector<std::pair<int, int>> forbiddenXY;
  std::string requiredRLE;
  std::pair<int, int> requiredXY;
  std::string antirequiredRLE;
  std::pair<int, int> antirequiredXY;
  std::string locusRLE;
  std::pair<int, int> locusXY;
  bool transparent;
  bool mustInclude;
  bool checkRecovery;
  bool sacrificial;
  unsigned period;
  bool isBlinker;
  unsigned maxUses = 0;  // at most this many per result (0 = no limit)
  std::vector<std::string> zoneRLE;            // this catalyst must lie inside these cells
  std::vector<std::pair<int, int>> zoneXY;

  explicit CatalystInput(std::string &line) {
    std::vector<std::string> elems = splitwhitespace(line);

    if (elems.size() < 6) {
      std::cout << "The line " << line << "is invalid" << std::endl;
      std::cout << "Format: cat <rle> <absense interval> <centerX> <centerY> "
                   "<symm Type | + / x *>"
                << std::endl;
      exit(1);
    }

    rle = elems[1];
    CheckRle(rle, "cat");
    maxDisappear = ParseIntOrDie(elems[2], "cat max-active", 0);
    centerX = ParseIntOrDie(elems[3], "cat dx");
    centerY = ParseIntOrDie(elems[4], "cat dy");
    symmType = elems[5].at(0);
    // Anything else used to be treated silently as '.', i.e. no extra orientations
    if (elems[5].size() != 1 || std::string(".|-\\/+@x*").find(symmType) == std::string::npos)
      Die("cat " + rle + ": symmetry \"" + elems[5] + "\" must be one of . | - \\ / + x *");

    transparent = false;
    mustInclude = false;
    checkRecovery = false;
    sacrificial = false;
    isBlinker = false;
    period = 1;

    unsigned argi = 6;

    // attributes that take "<rle> <x> <y>" or "<n>" need that many more tokens
    auto needTokens = [&](unsigned n) {
      if (argi + n >= elems.size())
        Die("cat " + rle + ": attribute \"" + elems[argi] + "\" is missing its arguments");
    };

    while (argi < elems.size()) {
      if (elems[argi] == "forbidden") {
        needTokens(3);
        CheckRle(elems[argi + 1], "cat forbidden");
        forbiddenRLE.push_back(elems[argi + 1]);
        forbiddenXY.emplace_back(
            ParseIntOrDie(elems[argi + 2], "cat x"), ParseIntOrDie(elems[argi + 3], "cat y"));

        argi += 4;
      } else if (elems[argi] == "required") {
        needTokens(3);
        CheckRle(elems[argi + 1], "cat required");
        requiredRLE = elems[argi + 1];
        requiredXY = std::make_pair(ParseIntOrDie(elems[argi + 2], "cat x"), ParseIntOrDie(elems[argi + 3], "cat y"));
        argi += 4;
      } else if (elems[argi] == "antirequired") {
        needTokens(3);
        CheckRle(elems[argi + 1], "cat antirequired");
        antirequiredRLE = elems[argi + 1];
        antirequiredXY = std::make_pair(ParseIntOrDie(elems[argi + 2], "cat x"), ParseIntOrDie(elems[argi + 3], "cat y"));
        argi += 4;
      } else if (elems[argi] == "locus") {
        needTokens(3);
        CheckRle(elems[argi + 1], "cat locus");
        locusRLE = elems[argi + 1];
        locusXY = std::make_pair(ParseIntOrDie(elems[argi + 2], "cat x"), ParseIntOrDie(elems[argi + 3], "cat y"));
        argi += 4;
      } else if (elems[argi] == "transparent") {
        transparent = true;
        argi += 1;
      } else if (elems[argi] == "mustinclude") {
        mustInclude = true;
        argi += 1;
      } else if (elems[argi] == "check-recovery") {
        checkRecovery = true;
        argi += 1;
      } else if (elems[argi] == "sacrificial") {
        sacrificial = true;
        argi += 1;
      } else if (elems[argi] == "isblinker") {
        isBlinker = true;
        argi += 1;
      } else if (elems[argi] == "max-uses") {
        needTokens(1);
        maxUses = ParseIntOrDie(elems[argi + 1], "cat max-uses", 1);
        argi += 2;
      } else if (elems[argi] == "zone") {
        needTokens(3);
        CheckRle(elems[argi + 1], "cat zone");
        zoneRLE.push_back(elems[argi + 1]);
        zoneXY.emplace_back(ParseIntOrDie(elems[argi + 2], "cat zone x"), ParseIntOrDie(elems[argi + 3], "cat zone y"));
        argi += 4;
      } else if (elems[argi] == "period") {
        needTokens(1);
        period = ParseIntOrDie(elems[argi + 1], "cat period", 1); // 0 would divide by zero later
        argi += 2;
      } else {
        std::cout << "Unknown catalyst attribute: " << elems[argi] << std::endl;
        exit(1);
      }
    }
  }

  void Print() const {
    std::cout << rle << " " << maxDisappear << " " << centerX << " " << centerY
              << " " << symmType << std::endl;
  }
};

inline std::vector<SymmetryTransform> SymmetryGroupFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {Identity};
  case StaticSymmetry::D2AcrossX:
    return {Identity, ReflectAcrossX};
    // vertical/horizontal here refer to box dimensions, NOT axis of reflection
  case StaticSymmetry::D2AcrossXEven:
    return {Identity, ReflectAcrossXEven};
  case StaticSymmetry::D2AcrossY:
    return {Identity, ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {Identity, ReflectAcrossYEven};
  case StaticSymmetry::D2diagodd:
    return {Identity, ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {Identity, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Identity, Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Identity, Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Identity, Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Identity, Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Identity, Rotate90, Rotate180OddBoth, Rotate270};
  case StaticSymmetry::C4even:
    return {Identity, Rotate90Even, Rotate180EvenBoth, Rotate270Even};
  case StaticSymmetry::D4:
    return {Identity, ReflectAcrossX, Rotate180OddBoth, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {Identity, ReflectAcrossXEven, Rotate180EvenBoth,
            ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {Identity, ReflectAcrossYEven, Rotate180EvenHorizontal,
            ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {Identity, ReflectAcrossXEven, Rotate180EvenVertical,
            ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {Identity, ReflectAcrossYeqX, Rotate180OddBoth,
            ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {Identity, ReflectAcrossYeqX, Rotate180EvenBoth,
            ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Identity,       ReflectAcrossX,         ReflectAcrossYeqX,
            ReflectAcrossY, ReflectAcrossYeqNegXP1, Rotate90,
            Rotate270,      Rotate180OddBoth};
  case StaticSymmetry::D8even:
    return {Identity,           ReflectAcrossXEven,   ReflectAcrossYeqX,
            ReflectAcrossYEven, ReflectAcrossYeqNegX, Rotate90Even,
            Rotate270Even,      Rotate180EvenBoth};
  }
}

inline std::vector<SymmetryTransform> SymmetryChainFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {};
  case StaticSymmetry::D2AcrossY:
    return {ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {ReflectAcrossYEven};
  case StaticSymmetry::D2AcrossX:
    return {ReflectAcrossX};
  case StaticSymmetry::D2AcrossXEven:
    return {ReflectAcrossXEven};
  case StaticSymmetry::D2diagodd:
    return {ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Rotate90, Rotate180OddBoth};
  case StaticSymmetry::C4even:
    return {Rotate90Even, Rotate180EvenBoth};
  case StaticSymmetry::D4:
    return {ReflectAcrossX, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {ReflectAcrossXEven, ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {ReflectAcrossYEven, ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {ReflectAcrossXEven, ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Rotate90, Rotate180OddBoth, ReflectAcrossYeqX};
  case StaticSymmetry::D8even:
    return {Rotate90Even, Rotate180EvenBoth, ReflectAcrossYeqX};
  }
}

// Cells (x, y) in storage coordinates 0..N-1 (the same cells the old fixed
// 64x64 RLE strings drew from the origin) where a catalyst may be placed.
LifeState FundamentalDomain(const StaticSymmetry sym) {
  auto build = [](bool (*inside)(int, int)) {
    LifeState r;
    for (int x = 0; x < N; x++)
      for (int y = 0; y < N; y++)
        if (inside(x, y))
          r.Set(x, y);
    r.RecalculateMinMax();
    return r;
  };
  const int H = N / 2;
  switch (sym) {
  case StaticSymmetry::C1:
    return build([](int, int) { return true; });
  case StaticSymmetry::D2AcrossY:
  case StaticSymmetry::D2AcrossYEven:
    return build([](int x, int) { return x < H; });
  case StaticSymmetry::D2AcrossX:
  case StaticSymmetry::D2AcrossXEven:
  case StaticSymmetry::C2:
  case StaticSymmetry::C2even:
  case StaticSymmetry::C2horizontaleven:
  case StaticSymmetry::C2verticaleven:
    return build([](int, int y) { return y < H; });
  case StaticSymmetry::D2diagodd:
    return build([](int x, int y) { return x + y <= N - 1; });
  case StaticSymmetry::D2negdiagodd:
    return build([](int x, int y) { return x <= y; });
  case StaticSymmetry::C4:
  case StaticSymmetry::C4even:
  case StaticSymmetry::D4:
  case StaticSymmetry::D4even:
  case StaticSymmetry::D4horizontaleven:
  case StaticSymmetry::D4verticaleven:
    return build([](int x, int y) { return x < H && y < H; });
  case StaticSymmetry::D4diag:
  case StaticSymmetry::D4diageven:
    return build([](int x, int y) { return x <= std::min(y, N - 1 - y); });
  case StaticSymmetry::D8:
  case StaticSymmetry::D8even:
    return build([](int x, int y) { return y < H && x <= y; });
  }
  return LifeState();
}

StaticSymmetry SymmetryFromString(const std::string &name) {
  std::string start = name.substr(0, 2);
  std::string rest = name.substr(2);
  if (start == "D2") {
    if (rest == "-" or rest == "vertical") {
      return StaticSymmetry::D2AcrossX;
    } else if (rest == "-even" or rest == "verticaleven") {
      return StaticSymmetry::D2AcrossXEven;
    } else if (rest == "|" or rest == "horizontal") {
      return StaticSymmetry::D2AcrossY;
    } else if (rest == "|even" or rest == "horizontaleven") {
      return StaticSymmetry::D2AcrossYEven;
    } else if (rest == "/" or rest == "/odd") {
      return StaticSymmetry::D2negdiagodd;
    } else if (rest == "\\" or rest == "\\odd") {
      return StaticSymmetry::D2diagodd;
    }
  } else if (start == "C2") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C2;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C2even;
    } else if (rest == "horizontaleven" or rest == "|even") {
      return StaticSymmetry::C2horizontaleven;
    } else if (rest == "verticaleven" or rest == "-even" or rest == "_2") {
      return StaticSymmetry::C2verticaleven;
    }
  } else if (start == "C4") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C4;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C4even;
    }
  } else if (start == "D4") {
    std::string evenOddInfo = rest.substr(1);
    if (rest[0] == '+' or (rest.size() > 1 and rest[1] == '+')) {
      if (evenOddInfo == "" or rest == "_+1") {
        return StaticSymmetry::D4;
      } else if (evenOddInfo == "even" or rest == "_+4") {
        return StaticSymmetry::D4even;
      } else if (evenOddInfo == "verticaleven" or evenOddInfo == "-even" or
                 rest == "_+2") {
        return StaticSymmetry::D4verticaleven;
      } else if (evenOddInfo == "horizontaleven" or evenOddInfo == "|even") {
        return StaticSymmetry::D4horizontaleven;
      }
    } else if (rest[0] == 'x' or (rest.size() > 1 and rest[1] == 'x')) {
      if (evenOddInfo == "" or rest == "_x1") {
        return StaticSymmetry::D4diag;
      } else if (evenOddInfo == "even" or rest == "_x4") {
        return StaticSymmetry::D4diageven;
      }
    }
  } else if (start == "D8") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::D8;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::D8even;
    }
  }
  return StaticSymmetry::C1;
}

std::vector<SymmetryTransform> CharToTransforms(char ch) {
  switch (ch) {
  case '.':
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  case '|':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossY);
  case '-':
    return SymmetryGroupFromEnum(StaticSymmetry::D2AcrossX);
  case '\\':
    return SymmetryGroupFromEnum(StaticSymmetry::D2diagodd);
  case '/':
    return SymmetryGroupFromEnum(StaticSymmetry::D2negdiagodd);
  case '+':
  case '@':
    return SymmetryGroupFromEnum(StaticSymmetry::C4);
  case 'x':
    return {Identity, Rotate90, ReflectAcrossX, ReflectAcrossYeqX};
  case '*':
    return SymmetryGroupFromEnum(StaticSymmetry::D8);
  default:
    return SymmetryGroupFromEnum(StaticSymmetry::C1);
  }
}

// Every keyword of the input format (anything else is warned about).
static const char *const kKeywords[] = {
    "rule", "include", "cat", "max-gen", "num-catalyst", "num-transparent", "stable-interval",
    "pat", "search-area", "exclude-area", "start-gen", "last-gen", "output", "full-report",
    "filter", "andfilter", "orfilter", "match", "fit-in-width-height", "max-category-size",
    "symmetry", "also-required", "stop-after-cats-destroyed", "max-junk", "keep-duplicates",
    "keep-interacting", "output-gap", "bait", "threads", "max-time", "max-results", "output-list",
    "verbose", "part", "unique-output", "lightcone", "lightcone-tight", "resume-from", "contact-filter", "split-level",
    "spaceship-cutoff", "torus", "fixed-cat", "catalyst-zone", "no-catalyst-zone", "keep-empty", "torus-edge", "sort-by",
    "layout", "legend"};

static bool IsKnownKeyword(const std::string &w) {
  for (const char *k : kKeywords)
    if (w == k)
      return true;
  return false;
}

static std::string Trim(const std::string &str) {
  size_t a = 0, b = str.size();
  while (a < b && isspace((unsigned char)str[a])) a++;
  while (b > a && isspace((unsigned char)str[b - 1])) b--;
  return str.substr(a, b - a);
}

// "x = 12, y = 5, rule = B3/S23" (an RLE header line)
static bool IsRleHeader(const std::string &t) {
  size_t i = 0;
  while (i < t.size() && isspace((unsigned char)t[i])) i++;
  if (i >= t.size() || t[i] != 'x') return false;
  i++;
  while (i < t.size() && isspace((unsigned char)t[i])) i++;
  return i < t.size() && t[i] == '=';
}

// Value of "key = value" in an RLE header (up to the next comma), or "".
static std::string HeaderField(const std::string &header, const std::string &key) {
  std::vector<std::string> parts;
  split(header, ',', parts);
  for (const std::string &part : parts) {
    size_t eq = part.find('=');
    if (eq != std::string::npos && Trim(part.substr(0, eq)) == key)
      return Trim(part.substr(eq + 1));
  }
  return "";
}

typedef std::map<std::pair<int, int>, int> CellStates; // (x, y) -> state, live cells only

// Decodes an RLE body (no header). b and . are state 0, o is 1, A..X are 1..24.
static CellStates DecodeRleBody(const std::string &body, const std::string &what) {
  CellStates cells;
  int x = 0, y = 0;
  long n = 0;
  bool haveCount = false;
  for (char ch : body) {
    if (isspace((unsigned char)ch))
      continue;
    if (isdigit((unsigned char)ch)) {
      n = n * 10 + (ch - '0');
      haveCount = true;
      if (n > 100000)
        Die(what + ": run length too large in RLE");
      continue;
    }
    int k = haveCount ? (int)n : 1;
    n = 0;
    haveCount = false;
    if (ch == '!')
      break;
    if (ch == '$') {
      y += k;
      x = 0;
      continue;
    }
    int v;
    if (ch == 'b' || ch == '.')
      v = 0;
    else if (ch == 'o')
      v = 1;
    else if (ch >= 'A' && ch <= 'X')
      v = ch - 'A' + 1;
    else
      Die(what + ": unsupported character '" + std::string(1, ch) + "' in RLE");
    if (v)
      for (int i = 0; i < k; i++)
        cells[{x + i, y}] = v;
    x += k;
  }
  return cells;
}

// Two-state RLE (o b $ !) of `cells`, drawn from the corner (x0, y0).
static std::string CellsToRle(const std::set<std::pair<int, int>> &cells, int x0, int y0) {
  if (cells.empty())
    return "!";
  int maxx = x0, maxy = y0;
  for (const auto &c : cells) {
    maxx = std::max(maxx, c.first);
    maxy = std::max(maxy, c.second);
  }
  std::string out;
  int blankRows = 0;
  for (int y = y0; y <= maxy; y++) {
    std::string row;
    int run = 0;
    bool cur = false;
    auto flush = [&]() {
      if (run > 0) {
        if (run > 1) row += std::to_string(run);
        row += cur ? 'o' : 'b';
      }
    };
    for (int x = x0; x <= maxx; x++) {
      bool v = cells.count({x, y}) > 0;
      if (run > 0 && v != cur) { flush(); run = 0; }
      cur = v;
      run++;
    }
    if (cur) flush(); // trailing dead cells are implied
    if (y > y0) blankRows++;
    if (!row.empty()) {
      if (blankRows > 0)
        out += (blankRows > 1 ? std::to_string(blankRows) : "") + "$";
      blankRows = 0;
      out += row;
    }
  }
  return out + "!";
}

// Same, cropped to the cells' bounding box; returns the box corner in x0, y0.
static std::string CellsToRleCropped(const std::set<std::pair<int, int>> &cells, int &x0, int &y0) {
  x0 = y0 = 0;
  if (cells.empty())
    return "!";
  x0 = cells.begin()->first;
  y0 = cells.begin()->second;
  for (const auto &c : cells) {
    x0 = std::min(x0, c.first);
    y0 = std::min(y0, c.second);
  }
  return CellsToRle(cells, x0, y0);
}

// 8-connected pieces of a cell set
static std::vector<std::set<std::pair<int, int>>> Components(const std::set<std::pair<int, int>> &cells) {
  std::vector<std::set<std::pair<int, int>>> out;
  std::set<std::pair<int, int>> left = cells;
  while (!left.empty()) {
    std::set<std::pair<int, int>> comp;
    std::vector<std::pair<int, int>> todo{*left.begin()};
    left.erase(left.begin());
    while (!todo.empty()) {
      auto c = todo.back();
      todo.pop_back();
      comp.insert(c);
      for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++) {
          auto it = left.find({c.first + dx, c.second + dy});
          if (it != left.end()) {
            todo.push_back(*it);
            left.erase(it);
          }
        }
    }
    out.push_back(comp);
  }
  return out;
}

// ---- layout files --------------------------------------------------------
// A layout is a multi-state RLE (LifeHistory states, so it can be drawn and
// edited in Golly) whose cell states give each cell a role; see README.md.
// It is turned into ordinary keyword lines (pat, bait, fixed-cat, ...).
struct LayoutRole {
  std::string kind; // pat zone exclude bait required fixed goal keep-empty ignore
  std::vector<std::string> args;
};

static void ApplyLegend(const std::string &spec, std::map<int, LayoutRole> &roles, const std::string &where) {
  std::vector<std::string> words = splitwhitespace(spec);
  for (size_t i = 0; i < words.size(); i++) {
    if (i == 0 && words[i] == "legend")
      continue;
    size_t eq = words[i].find('=');
    if (eq == std::string::npos)
      Die(where + ": legend entry \"" + words[i] + "\" must look like 4=goal");
    int state = ParseIntOrDie(words[i].substr(0, eq), where + ": legend state", 1);
    if (state > 24)
      Die(where + ": legend state " + std::to_string(state) + " is above 24");
    std::vector<std::string> parts;
    split(words[i].substr(eq + 1), ':', parts);
    if (parts.empty())
      Die(where + ": legend entry \"" + words[i] + "\" has no role");
    std::string kind = parts[0];
    if (kind == "catalyst-zone") kind = "zone";
    if (kind == "no-catalyst-zone") kind = "exclude";
    if (kind == "also-required") kind = "required";
    if (kind == "fixed-cat") kind = "fixed";
    static const char *const kinds[] = {"pat", "zone", "exclude", "bait", "required", "fixed", "goal", "keep-empty", "ignore"};
    if (std::find_if(std::begin(kinds), std::end(kinds), [&](const char *k) { return kind == k; }) == std::end(kinds))
      Die(where + ": unknown legend role \"" + kind + "\" (use pat, zone, exclude, bait, required, fixed, goal, keep-empty, ignore)");
    roles[state] = {kind, std::vector<std::string>(parts.begin() + 1, parts.end())};
  }
}

// Keyword lines for the layout in `fname`. Its top-left corner goes to
// (ox, oy) if given, else to the "#CXRLE Pos=" Golly writes, else the layout
// is centred on (0, 0). `outerLegend`: legend lines of the including file.
static std::vector<std::string> LayoutToLines(const std::string &fname, bool originGiven, int ox, int oy,
                                              const std::vector<std::string> &outerLegend) {
  std::ifstream in(fname.c_str());
  if (!in.good())
    Die("Could not open layout file: " + fname);

  std::vector<std::string> keywordLines, fileLegend;
  std::string header, body, line;
  bool havePos = false;
  int px = 0, py = 0;
  while (std::getline(in, line)) {
    std::string t = Trim(line);
    if (t.empty())
      continue;
    if (t.compare(0, 6, "#CXRLE") == 0) {
      size_t pos = t.find("Pos=");
      if (pos != std::string::npos && sscanf(t.c_str() + pos + 4, "%d,%d", &px, &py) == 2)
        havePos = true;
    } else if (t.size() >= 2 && t[0] == '#' && (t[1] == 'C' || t[1] == 'c')) {
      // "#C <keyword> ..." lines are part of the search; other comments are ignored
      std::string rest = Trim(t.substr(2));
      std::vector<std::string> words = splitwhitespace(rest);
      if (words.empty())
        continue;
      if (words[0] == "legend")
        fileLegend.push_back(rest);
      else if (IsKnownKeyword(words[0]))
        keywordLines.push_back(rest);
    } else if (t[0] == '#') {
      continue;
    } else if (header.empty() && body.empty() && IsRleHeader(t)) {
      header = t;
    } else if (body.find('!') != std::string::npos) {
      keywordLines.push_back(t); // after the RLE: plain keyword lines
    } else {
      body += t;
    }
  }

  std::map<int, LayoutRole> roles;
  ApplyLegend("1=pat 2=zone 3=bait 4=goal 5=fixed 6=keep-empty", roles, fname);
  for (const std::string &l : fileLegend) ApplyLegend(l, roles, fname);
  for (const std::string &l : outerLegend) ApplyLegend(l, roles, fname);

  CellStates cells = DecodeRleBody(body, "layout " + fname);
  int w = 0, h = 0;
  for (const auto &c : cells) {
    w = std::max(w, c.first.first + 1);
    h = std::max(h, c.first.second + 1);
  }
  if (!header.empty()) {
    std::string hw = HeaderField(header, "x"), hh = HeaderField(header, "y");
    if (!hw.empty()) w = std::max(w, atoi(hw.c_str()));
    if (!hh.empty()) h = std::max(h, atoi(hh.c_str()));
  }
  if (w > N || h > N)
    Die("layout " + fname + ": " + std::to_string(w) + "x" + std::to_string(h) + " does not fit on the " +
        std::to_string(N) + "x" + std::to_string(N) + " torus" + (N < 128 ? " (add \"torus 128\")" : ""));
  if (!originGiven) {
    if (havePos) { ox = px; oy = py; }
    else { ox = -(w / 2); oy = -(h / 2); }
  }

  std::vector<std::string> out;
  out.push_back("# layout " + fname + ": its top-left cell is at (" + std::to_string(ox) + "," + std::to_string(oy) + ")");
  std::string rule = HeaderField(header, "rule");
  for (const char *suffix : {"History", "history", "Super", "super"}) {
    size_t L = strlen(suffix);
    if (rule.size() > L && rule.compare(rule.size() - L, L, suffix) == 0)
      rule = rule.substr(0, rule.size() - L);
  }
  if (!rule.empty())
    out.push_back("rule " + rule);

  std::map<std::string, std::set<std::pair<int, int>>> byKind;
  std::map<int, std::set<std::pair<int, int>>> byState;
  for (const auto &c : cells) {
    auto r = roles.find(c.second);
    if (r == roles.end())
      Die("layout " + fname + ": cells in state " + std::to_string(c.second) +
          " have no meaning; give it one with a legend line, e.g. \"#C legend " + std::to_string(c.second) + "=ignore\"");
    byState[c.second].insert({c.first.first + ox, c.first.second + oy});
  }
  auto arg = [](const LayoutRole &r, size_t i) { return i < r.args.size() ? r.args[i] : std::string(); };
  for (const auto &st : byState) {
    const LayoutRole &r = roles[st.first];
    const auto &set = st.second;
    int x0, y0;
    if (r.kind == "ignore")
      continue;
    if (r.kind == "pat" || r.kind == "required" || r.kind == "zone" || r.kind == "exclude" || r.kind == "keep-empty") {
      byKind[r.kind].insert(set.begin(), set.end());
      continue;
    }
    for (const auto &comp : Components(set)) {
      std::string rle = CellsToRleCropped(comp, x0, y0);
      std::string pos = " " + std::to_string(x0) + " " + std::to_string(y0);
      if (r.kind == "bait") {
        std::string l = "bait " + rle + pos;
        if (!arg(r, 0).empty()) l += " survive " + arg(r, 0);
        if (!arg(r, 1).empty()) l += " absent " + arg(r, 1);
        out.push_back(l);
      } else if (r.kind == "fixed") {
        std::string l = "fixed-cat " + rle + pos + " " + (arg(r, 0).empty() ? "20" : arg(r, 0));
        for (size_t i = 1; i < r.args.size(); i++) l += " " + r.args[i];
        out.push_back(l);
      } else if (r.kind == "goal") {
        out.push_back("andfilter " + (arg(r, 0).empty() ? std::string("1-max") : arg(r, 0)) + " " + rle + pos);
      }
    }
  }
  static const std::pair<const char *, const char *> unions[] = {
      {"pat", "pat"}, {"required", "also-required"}, {"zone", "catalyst-zone"},
      {"exclude", "no-catalyst-zone"}, {"keep-empty", "keep-empty"}};
  for (const auto &u : unions) {
    auto it = byKind.find(u.first);
    if (it == byKind.end())
      continue;
    int x0, y0;
    std::string rle = CellsToRleCropped(it->second, x0, y0);
    out.push_back(std::string(u.second) + " " + rle + " " + std::to_string(x0) + " " + std::to_string(y0));
  }
  out.insert(out.end(), keywordLines.begin(), keywordLines.end());
  return out;
}

// True if the file is a layout (an RLE) rather than a keyword file.
static bool IsLayoutFile(const std::string &fname) {
  if (fname.size() > 4 && (fname.compare(fname.size() - 4, 4, ".rle") == 0 || fname.compare(fname.size() - 4, 4, ".RLE") == 0))
    return true;
  std::ifstream in(fname.c_str());
  std::string line;
  while (std::getline(in, line)) {
    std::string t = Trim(line);
    if (t.empty() || t[0] == '#')
      continue;
    return IsRleHeader(t);
  }
  return false;
}

static std::vector<std::string> ReadInputLines(const std::string &fname, int depth = 0);

// Turns raw lines into keyword lines: strips "  # comments" at line ends,
// reads "<<" RLE blocks, and expands include and layout lines.
static std::vector<std::string> ProcessInputLines(const std::vector<std::string> &raw, const std::string &dir, int depth) {
  auto resolve = [&](const std::string &path) {
    if (!dir.empty() && std::ifstream((dir + path).c_str()).good())
      return dir + path;
    return path;
  };
  // a '#' after a space starts a comment (on keyword lines)
  auto stripComment = [](std::string line) {
    if (!line.empty() && line.back() == '\r')
      line.pop_back();
    std::string t = Trim(line);
    if (!t.empty() && isalpha((unsigned char)t[0]))
      for (size_t i = 1; i < line.size(); i++)
        if (line[i] == '#' && isspace((unsigned char)line[i - 1]))
          return line.substr(0, i);
    return line;
  };
  std::vector<std::string> legends;
  for (const std::string &l : raw) {
    std::string c = stripComment(l);
    std::vector<std::string> e = splitwhitespace(c);
    if (!e.empty() && e[0] == "legend")
      legends.push_back(Trim(c));
  }

  std::vector<std::string> out;
  for (size_t li = 0; li < raw.size(); li++) {
    std::string line = stripComment(raw[li]);
    std::vector<std::string> elems = splitwhitespace(line);

    // "<<" where an RLE goes: the RLE follows on the next lines, up to its "!"
    // (as copied from Golly: header and #-lines are skipped)
    auto here = std::find(elems.begin(), elems.end(), "<<");
    if (here != elems.end()) {
      std::string body;
      bool closed = false;
      while (++li < raw.size()) {
        std::string r = Trim(raw[li]);
        if (r.empty() || r[0] == '#' || IsRleHeader(r))
          continue;
        body += r;
        if (r.find('!') != std::string::npos) { closed = true; break; }
      }
      if (!closed)
        Die("The RLE after \"" + Trim(line) + "\" has no closing '!'");
      CellStates cells = DecodeRleBody(body, elems[0]);
      std::set<std::pair<int, int>> live;
      for (const auto &c : cells)
        if (c.second % 2 == 1) // LifeHistory: odd states are alive
          live.insert(c.first);
      *here = CellsToRle(live, 0, 0); // from the RLE's own corner, so offsets keep their meaning
      line.clear();
      for (const std::string &e : elems)
        line += (line.empty() ? "" : " ") + e;
    }

    if (elems.size() >= 2 && elems[0] == "include") {
      std::string path = Trim(Trim(line).substr(7));
      std::vector<std::string> sub = ReadInputLines(resolve(path), depth + 1);
      out.insert(out.end(), sub.begin(), sub.end());
    } else if (elems.size() >= 2 && elems[0] == "layout") {
      if (elems.size() != 2 && elems.size() != 4)
        Die("Format: layout <file.rle> [x y]");
      bool given = elems.size() == 4;
      int x = given ? ParseIntOrDie(elems[2], "layout x") : 0, y = given ? ParseIntOrDie(elems[3], "layout y") : 0;
      std::vector<std::string> sub = LayoutToLines(resolve(elems[1]), given, x, y, legends);
      sub = ProcessInputLines(sub, dir, depth + 1);
      out.insert(out.end(), sub.begin(), sub.end());
    } else if (!elems.empty() && elems[0] == "legend") {
      continue; // used by the layout lines above
    } else {
      out.push_back(line);
    }
  }
  return out;
}

// Reads an input file: a keyword file (.in) or a layout (.rle). Expands
// "include <file>" (relative paths are tried next to the including file first,
// then in the current directory), "layout <file>", "<<" RLE blocks and end-of-line
// comments. Also strips Windows line endings.
static std::vector<std::string> ReadInputLines(const std::string &fname, int depth) {
  if (depth > 8)
    Die("include: nested too deeply (is there a loop?) at " + fname);

  std::string dir;
  size_t slash = fname.find_last_of("/\\");
  if (slash != std::string::npos)
    dir = fname.substr(0, slash + 1);

  if (IsLayoutFile(fname))
    return ProcessInputLines(LayoutToLines(fname, false, 0, 0, {}), dir, depth);

  std::ifstream in(fname.c_str());
  if (!in.good())
    Die("Could not open file: " + fname);
  std::vector<std::string> raw;
  std::string line;
  while (std::getline(in, line))
    raw.push_back(line);
  return ProcessInputLines(raw, dir, depth);
}

// Scans the input for a "rule" line and applies it. This is done in a
// separate first pass so the rule is set before any pattern/catalyst is
// stepped, regardless of where the line appears in the file.
// Default (no rule line): B3/S23 (Conway's Game of Life).
static void ApplyRuleFromFile(const std::vector<std::string> &lines) {
  std::string ruleStr;
  for (const std::string &line : lines) {
    std::vector<std::string> elems = splitwhitespace(line);
    if (elems.size() >= 2 && elems[0] == "rule") {
      ruleStr.clear();
      for (unsigned k = 1; k < elems.size(); k++)
        ruleStr += elems[k];
    }
  }
  if (ruleStr.empty())
    ruleStr = "B3/S23";
  std::string err;
  if (!SetRule(ruleStr, err)) {
    std::cout << "Invalid rule \"" << ruleStr << "\": " << err << std::endl;
    exit(1);
  }
  std::cout << "Rule: " << RuleDescription() << std::endl;
}

void ReadParams(const std::string& fname, std::vector<CatalystInput> &catalysts,
                SearchParams &params) {
  std::vector<std::string> lines = ReadInputLines(fname);

  ApplyRuleFromFile(lines);

  std::string Cat = "cat";
  std::string maxGen = "max-gen";
  std::string startGen = "start-gen";
  std::string lastGen = "last-gen";

  std::string numCat = "num-catalyst";
  std::string numTransp = "num-transparent";
  std::string stable = "stable-interval";
  std::string area = "search-area";
  std::string pat = "pat";
  std::string outputFile = "output";
  std::string filter = "filter";
  std::string andfilter = "andfilter";
  std::string orfilter = "orfilter";
  std::string matchfilter = "match";
  std::string maxWH = "fit-in-width-height";
  std::string maxCatSize = "max-category-size";
  std::string fullReport = "full-report";

  std::string symmetry = "symmetry";
  std::string alsoRequired = "also-required";
  std::string stopAfterCatsDestroyed = "stop-after-cats-destroyed";
  std::string maxJunk = "max-junk";
  std::string excludeArea = "exclude-area";
  std::string baitKw = "bait";
  std::string keepDup = "keep-duplicates";
  std::string outGap = "output-gap";

  std::string line;

  bool badSymmetry = false;
  bool hasLastGen = false;

  for (size_t lineNo = 0; lineNo < lines.size(); lineNo++) {
    line = lines[lineNo];
    std::vector<std::string> elems = splitwhitespace(line);

    // Lines starting with anything but a letter (#, =, -, ...) are comments.
    if (elems.empty() || !isalpha((unsigned char)elems[0][0]))
      continue;
    bool isKnown = IsKnownKeyword(elems[0]);
    if (!isKnown) {
      std::cout << "Warning: unknown keyword \"" << elems[0] << "\" (line ignored): " << line << std::endl;
      continue;
    }
    if (elems.size() < 2) {
      std::cout << "Warning: \"" << elems[0] << "\" needs a value (line ignored)" << std::endl;
      continue;
    }

    if (elems[0] == "include")
      continue; // already expanded by ReadInputLines

    if (elems[0] == "rule")
      continue; // already applied by ApplyRuleFromFile before anything is stepped

    if (elems[0] == Cat)
      catalysts.emplace_back(line);

    if (elems[0] == maxGen)
      params.maxGen = ParseIntOrDie(elems[1], "max-gen", 1);

    if (elems[0] == numCat) {
      // The search keeps per-catalyst state in fixed arrays of MAX_CATALYSTS entries.
      // "num-catalyst 1-3": every result with 1 to 3 catalysts.
      std::vector<std::string> range;
      split(elems[1], '-', range);
      if (range.size() == 2) {
        params.minCatalysts = ParseIntOrDie(range[0], "num-catalyst (fewest)", 1);
        params.numCatalysts = ParseIntOrDie(range[1], "num-catalyst (most)", 1);
        if (params.minCatalysts > params.numCatalysts)
          Die("num-catalyst " + elems[1] + ": the first number must not be larger than the second");
      } else {
        params.numCatalysts = ParseIntOrDie(elems[1], "num-catalyst", 1);
        params.minCatalysts = params.numCatalysts;
      }
      if (params.numCatalysts > (unsigned)MAX_SEARCHED)
        Die("num-catalyst " + elems[1] + " is not supported (maximum " + std::to_string(MAX_SEARCHED) + ")");
    }

    if (elems[0] == numTransp)
      params.numTransparent = ParseIntOrDie(elems[1], "num-transparent", 0);

    if (elems[0] == stable)
      params.stableInterval = ParseIntOrDie(elems[1], "stable-interval", 0);

    if (elems[0] == pat) {
      if (params.patSeen)
        Die("pat: the active pattern is given twice (a layout file also sets it)");
      params.patSeen = true;
      CheckRle(elems[1], "pat");
      params.pat = elems[1];

      if (elems.size() > 3) {
        params.xPat = ParseIntOrDie(elems[2], "pat x");
        params.yPat = ParseIntOrDie(elems[3], "pat y");
      }
    }

    if (elems[0] == area || elems[0] == excludeArea) {
      if (elems.size() < 5) {
        std::cout << "Format: " << elems[0] << " <x> <y> <width> <height>" << std::endl;
        exit(1);
      }
      std::array<int, 4> r = {ParseIntOrDie(elems[1], elems[0] + " x"), ParseIntOrDie(elems[2], elems[0] + " y"),
                              ParseIntOrDie(elems[3], elems[0] + " width"), ParseIntOrDie(elems[4], elems[0] + " height")};
      if (r[2] <= 0 || r[3] <= 0) {
        std::cout << elems[0] << ": width and height must be positive" << std::endl;
        exit(1);
      }
      if (elems[0] == area) {
        params.searchAreas.push_back(r); // several lines are combined (union)
        // keep the single-rectangle fields in sync with the first one
        if (params.searchAreas.size() == 1)
          for (int k = 0; k < 4; k++)
            params.searchArea[k] = r[k];
      } else {
        params.excludeAreas.push_back(r);
      }
    }

    if (elems[0] == startGen)
      params.startGen = ParseIntOrDie(elems[1], "start-gen", 0);

    if (elems[0] == lastGen) {
      params.lastGen = ParseIntOrDie(elems[1], "last-gen", 0);
      hasLastGen = true;
    }

    if (elems[0] == outputFile) {
      params.outputFile = elems[1];

      for (unsigned i = 2; i < elems.size(); i++) {
        params.outputFile.append(" ");
        params.outputFile.append(elems[i]);
      }
    }

    if (elems[0] == fullReport) {
      params.fullReportFile = elems[1];

      for (unsigned i = 2; i < elems.size(); i++) {
        params.fullReportFile.append(" ");
        params.fullReportFile.append(elems[i]);
      }
    }

    if (elems[0] == filter || elems[0] == orfilter || elems[0] == andfilter || elems[0] == matchfilter) {
      const bool isMatch = elems[0] == matchfilter;
      if (elems.size() < (isMatch ? 3u : 5u))
        Die("Format: " + elems[0] + (isMatch ? " <gen or range> <rle>" : " <gen or range> <rle> <dx> <dy>"));
      CheckRle(elems[2], elems[0]);

      std::vector<std::string> rangeElems;
      split(elems[1], '-', rangeElems);
      if (rangeElems.empty() || rangeElems.size() > 2)
        Die(elems[0] + ": generation \"" + elems[1] + "\" must be like 50 or 1-100");

      if (rangeElems.size() == 1) {
        params.filterGen.push_back(ParseIntOrDie(elems[1], elems[0] + " generation", 0));
        params.filterGenRange.emplace_back(-1, -1);
      } else {
        unsigned minGen = ParseIntOrDie(rangeElems[0], elems[0] + " generation range", 0);
        // "max" = the last generation (max-gen - 1), filled in after all lines are read
        unsigned maxGen = rangeElems[1] == "max" ? kFilterMaxGen : ParseIntOrDie(rangeElems[1], elems[0] + " generation range", 0);

        params.filterGen.push_back(-1);
        params.filterGenRange.emplace_back(minGen, maxGen);
      }

      params.targetFilter.push_back(elems[2]);

      if (elems[0] == andfilter) {
        params.filterdx.push_back(ParseIntOrDie(elems[3], elems[0] + " dx"));
        params.filterdy.push_back(ParseIntOrDie(elems[4], elems[0] + " dy"));
        params.filterType.push_back(ANDFILTER);
      } else if (elems[0] == orfilter) {
        params.filterdx.push_back(ParseIntOrDie(elems[3], elems[0] + " dx"));
        params.filterdy.push_back(ParseIntOrDie(elems[4], elems[0] + " dy"));
        params.filterType.push_back(ORFILTER);
      } else if (elems[0] == matchfilter) {
        params.filterdx.push_back(0);
        params.filterdy.push_back(0);
        params.filterType.push_back(MATCHFILTER);
      } else {
        // Shouldn't happen
        params.filterdx.push_back(ParseIntOrDie(elems[3], elems[0] + " dx"));
        params.filterdy.push_back(ParseIntOrDie(elems[4], elems[0] + " dy"));
        params.filterType.push_back(ANDFILTER);
      }
    }

    if (elems[0] == maxWH) {
      if (elems.size() < 3)
        Die("Format: fit-in-width-height <w> <h>");
      params.maxW = ParseIntOrDie(elems[1], "fit-in-width-height w", 1);
      params.maxH = ParseIntOrDie(elems[2], "fit-in-width-height h", 1);
    }

    if (elems[0] == maxCatSize) {
      params.maxCatSize = atoi(elems[1].c_str());
    }

    if (elems[0] == symmetry) {
      // Old spellings (the original aliases "horizontal" and "rotate180evenboth"
      // were mapped to names the parser did not know, i.e. silently no symmetry).
      std::string name = elems[1];
      if (name == "horizontal") name = "D2|";
      else if (name == "horizontaleven") name = "D2|even";
      else if (name == "diagonal") name = "D2/";
      else if (name == "rotate180") name = "C2";
      else if (name == "rotate180evenx") name = "C2|even";
      else if (name == "rotate180evenboth") name = "C2even";
      else if (name == "none") name = "C1";

      StaticSymmetry sym = StaticSymmetry::C1;
      bool ok = true;
      try {
        sym = SymmetryFromString(name);
      } catch (const std::exception &) { // e.g. "D4" made the parser throw
        ok = false;
      }
      // SymmetryFromString returns C1 for anything it does not recognise
      if (ok && sym == StaticSymmetry::C1 && name != "C1")
        ok = false;
      if (!ok)
        Die("Unknown symmetry \"" + elems[1] + "\". Use one of: C1, D2|, D2-, D2\\, D2/, D2|even, D2-even, "
            "C2, C2even, C2|even, C2-even, C4, C4even, D4+, D4+even, D4+|even, D4+-even, D4x, D4xeven, D8, D8even");

      params.symmetry = sym;
      params.symmetryName = name;
      params.symmetryChain = SymmetryChainFromEnum(params.symmetry);
    }
    if (elems[0] == alsoRequired) {
      if (elems.size() < 4)
        Die("Format: also-required <rle> <x> <y>");
      CheckRle(elems[1], "also-required");
      params.alsoRequired = elems[1].c_str();
      params.alsoRequiredXY = std::make_pair(ParseIntOrDie(elems[2], "also-required x"), ParseIntOrDie(elems[3], "also-required y"));
    }
    if (elems[0] == stopAfterCatsDestroyed){
      params.stopAfterCatsDestroyed = atoi(elems[1].c_str());
    }
    if (elems[0] == maxJunk){
      params.maxJunk = atoi(elems[1].c_str());
    }
    if (elems[0] == "keep-interacting") {
      std::string v = elems[1];
      params.keepInteracting = !(v == "no" || v == "false" || v == "0" || v == "off");
    }
    if (elems[0] == "verbose") {
      std::string v = elems[1];
      params.verbose = !(v == "no" || v == "false" || v == "0" || v == "off");
    }
    if (elems[0] == "resume-from")
      params.resumeFrom = std::stoul(elems[1]);
    if (elems[0] == "split-level")
      params.splitLevel = std::stoul(elems[1]);
    if (elems[0] == "contact-filter")
      params.soloRecovery = !(elems[1] == "no" || elems[1] == "false" || elems[1] == "0" || elems[1] == "off");
    if (elems[0] == "lightcone-tight")
      params.lightconeTight = !(elems[1] == "no" || elems[1] == "false" || elems[1] == "0" || elems[1] == "off");
    if (elems[0] == "lightcone")
      params.lightcone = !(elems[1] == "no" || elems[1] == "false" || elems[1] == "0" || elems[1] == "off");
    if (elems[0] == "unique-output")
      params.uniqueFile = elems[1];
    if (elems[0] == "fixed-cat" || elems[0] == "catalyst-zone" || elems[0] == "no-catalyst-zone" ||
        elems[0] == "keep-empty") {
      const bool fixed = elems[0] == "fixed-cat";
      if (elems.size() < (fixed ? 5u : 4u))
        Die("Format: " + elems[0] + (fixed ? " <rle> <x> <y> <max-active> [cat options]" : " <rle> <x> <y>"));
      CheckRle(elems[1], elems[0]);
      SearchParams::PlacedRle r;
      r.rle = elems[1];
      r.x = ParseIntOrDie(elems[2], elems[0] + " x");
      r.y = ParseIntOrDie(elems[3], elems[0] + " y");
      if (fixed) {
        ParseIntOrDie(elems[4], "fixed-cat max-active", 0);
        for (size_t k = 4; k < elems.size(); k++)
          r.rest += (k > 4 ? " " : "") + elems[k];
        params.fixedCats.push_back(r);
      } else {
        (elems[0] == "catalyst-zone" ? params.catalystZones : elems[0] == "keep-empty" ? params.keepEmpty
                                                                                         : params.noCatalystZones).push_back(r);
      }
    }
    if (elems[0] == "sort-by") {
      std::vector<std::string> keys;
      split(elems[1], ',', keys);
      for (const std::string &k : keys)
        if (k != "found" && k != "recovery" && k != "cells" && k != "catalysts" && k != "junk" && k != "frequency")
          Die("sort-by: unknown key \"" + k + "\" (use found, recovery, cells, catalysts, junk, frequency; several separated by commas)");
      params.sortBy = elems[1];
    }
    if (elems[0] == "torus-edge") {
      params.torusEdge = elems[1];
      if (elems[1] != "flag" && elems[1] != "drop" && elems[1] != "off")
        Die("torus-edge: use flag, drop or off");
    }
    if (elems[0] == "spaceship-cutoff") {
      const std::string &v = elems[1];
      params.spaceshipCutoff = (v == "no" || v == "off" || v == "false") ? 0 : ParseIntOrDie(v, "spaceship-cutoff", 0);
    }
    if (elems[0] == "torus") {
      params.torus = ParseIntOrDie(elems[1], "torus", 1);
      if (params.torus != 64 && params.torus != 128)
        Die("torus: " + elems[1] + " is not supported (use 64 or 128)");
    }
    if (elems[0] == "part")
      ParsePart(elems[1], params.partIndex, params.numParts);
    if (elems[0] == "threads")
      params.numThreads = ParseIntOrDie(elems[1], "threads", 0);
    if (elems[0] == "max-time")
      params.maxTimeSec = ParseDuration(elems[1], "max-time");
    if (elems[0] == "max-results")
      params.maxResults = ParseIntOrDie(elems[1], "max-results", 1);
    if (elems[0] == "output-list") {
      params.outputList = elems[1];
      for (unsigned i = 2; i < elems.size(); i++)
        params.outputList += " " + elems[i];
    }
    if (elems[0] == baitKw) {
      if (elems.size() < 4) {
        std::cout << "Format: bait <rle> <x> <y> [survive <gens>] [absent <gens>]" << std::endl;
        exit(1);
      }
      SearchParams::BaitInput b;
      CheckRle(elems[1], "bait");
      b.rle = elems[1];
      b.x = ParseIntOrDie(elems[2], "bait x");
      b.y = ParseIntOrDie(elems[3], "bait y");
      for (unsigned k = 4; k < elems.size(); k++) {
        if (elems[k] == "survive" && k + 1 < elems.size())
          b.survive = ParseIntOrDie(elems[++k], "bait survive", 0);
        else if (elems[k] == "absent" && k + 1 < elems.size())
          b.absent = ParseIntOrDie(elems[++k], "bait absent", 0);
        else {
          std::cout << "Unknown bait attribute \"" << elems[k] << "\"" << std::endl;
          exit(1);
        }
      }
      if (params.baits.size() >= (size_t)MAX_BAITS) {
        std::cout << "Too many bait lines (max " << MAX_BAITS << ")" << std::endl;
        exit(1);
      }
      params.baits.push_back(b);
    }
    if (elems[0] == keepDup) {
      std::string v = elems[1];
      params.keepDuplicates = !(v == "no" || v == "false" || v == "0" || v == "off");
    }
    if (elems[0] == outGap) {
      params.outputGap = std::max(0, atoi(elems[1].c_str()));
    }
  }

  if(!hasLastGen)
    params.lastGen = params.maxGen - 1;
  if (params.minCatalysts == 0)
    params.minCatalysts = params.numCatalysts;
  for (auto &r : params.filterGenRange)
    if (r.second == (int)kFilterMaxGen)
      r.second = (int)params.maxGen - 1;

  if (params.pat.length() == 0) {
    std::cout << "Did not read any pattern!" << std::endl;
    exit(1);
  }
  if (catalysts.empty()) {
    std::cout << "Did not read any catalysts!" << std::endl;
    exit(1);
  }
  if (badSymmetry) {
    std::cout << "Couldn\'t parse symmetry option" << std::endl;
    exit(1);
  }
}

// Per-thread work partition: with N threads, thread t explores the top-level
// placements whose running index is congruent to t modulo N.
static thread_local unsigned tlsPart = 0;
static thread_local unsigned tlsParts = 1;
static thread_local unsigned long tlsRootCounter = 0;
// Pause/resume: top-level placement each thread is working on (all lower ones are finished).
static std::unique_ptr<std::atomic<unsigned long>[]> gCurRoot;
static unsigned gNumWorkers = 0;
// Dynamic load balancing: threads claim top-level placements in order from a shared ticket.
static std::atomic<unsigned long> gNextTicket{0};
static thread_local unsigned long tlsTicket = ~0UL, tlsQual = 0;

// Set by Ctrl-C; the search then stops cleanly and still saves its results.
static volatile std::sig_atomic_t gInterrupted = 0;
static void OnInterrupt(int) {
  if (gInterrupted)
    std::_Exit(130); // second Ctrl-C: quit immediately
  gInterrupted = 1;
}

// RLE (with "!") of the live cells of `st` inside their bounding box; also
// returns the box (top-left x0,y0 in torus coordinates, size w x h).
static std::string BoundingBoxRle(const LifeState &st, int &x0, int &y0, int &w, int &h) {
  int minx = N, miny = N, maxx = -N, maxy = -N;
  for (int y = -N / 2; y < N / 2; y++)
    for (int x = -N / 2; x < N / 2; x++)
      if (st.GetCell(x, y) == 1) {
        minx = std::min(minx, x); maxx = std::max(maxx, x);
        miny = std::min(miny, y); maxy = std::max(maxy, y);
      }
  if (maxx < minx) {
    x0 = y0 = w = h = 0;
    return "!";
  }
  x0 = minx; y0 = miny; w = maxx - minx + 1; h = maxy - miny + 1;
  std::string out;
  for (int y = miny; y <= maxy; y++) {
    std::string row;
    int run = 0;
    bool cur = false;
    auto flush = [&]() {
      if (run > 0) {
        if (run > 1) row += std::to_string(run);
        row += cur ? 'o' : 'b';
      }
    };
    for (int x = minx; x <= maxx; x++) {
      bool v = st.GetCell(x, y) == 1;
      if (run > 0 && v != cur) { flush(); run = 0; }
      cur = v; run++;
    }
    if (cur) flush(); // trailing dead cells are implied
    if (y > miny) out += '$';
    out += row;
  }
  return out + "!";
}

// A catalyst's recovery condition: satisfied if ANY alternative LifeTarget
// matches. Normal catalysts have exactly one alternative (their own phase);
// "isblinker" catalysts have one per phase of the oscillator, so they count
// as recovered in any phase.
struct TargetSet {
  std::vector<LifeTarget> alts;
  std::vector<row_t> altCols; // columns where alts[i] has wanted or unwanted cells

  void Add(const LifeTarget &t) {
    alts.push_back(t);
    row_t cols = 0;
    for (int i = 0; i < N; i++)
      if ((t.wanted.state[i] | t.unwanted.state[i]) != 0)
        cols |= (row_t)1 << i;
    altCols.push_back(cols);
  }

  void Clear() {
    alts.clear();
    altCols.clear();
  }

  // Same result as state.Contains(alt) for some alt, but only the columns the
  // target occupies are looked at (a catalyst spans a handful of the 64).
  bool MatchedBy(const LifeState &state) const {
    for (size_t k = 0; k < alts.size(); k++) {
      const LifeTarget &t = alts[k];
      bool ok = true;
      for (row_t bits = altCols[k]; bits; bits &= bits - 1) {
        int i = RowCtz(bits);
        row_t have = state.state[i];
        if ((t.wanted.state[i] & ~have) != 0 || (t.unwanted.state[i] & have) != 0) {
          ok = false;
          break;
        }
      }
      if (ok)
        return true;
    }
    return false;
  }
};

class CatalystData {
public:
  LifeState state;
  std::vector<LifeState> phases;
  LifeTarget target;
  LifeState reactionMask;
  std::vector<LifeState> phaseReactionMask;
  unsigned maxDisappear;
  std::vector<LifeTarget> forbidden;
  LifeState required;
  LifeState antirequired;
  bool hasLocus;
  LifeState locus;
  LifeState locusReactionMask;
  LifeState locusAvoidMask;
  bool transparent;
  bool mustInclude;
  bool checkRecovery;
  bool sacrificial;
  unsigned period;
  bool isBlinker;
  std::string name;      // the RLE from its cat line
  int inputIndex = -1;   // which cat line (0-based) it came from
  int orient = 0;        // which orientation of that line
  bool fixed = false;    // a fixed-cat: already placed, never searched for
  unsigned maxUses = 0;  // cat option max-uses (0 = no limit)
  LifeState zone;        // cat option zone: must lie inside these cells (empty = anywhere)

  static std::vector<CatalystData> FromInput(CatalystInput &input);
};

std::vector<CatalystData> CatalystData::FromInput(CatalystInput &input) {
  std::vector<SymmetryTransform> trans = CharToTransforms(input.symmType);

  const char *rle = input.rle.c_str();

  std::vector<CatalystData> results;
  static int nextInputIndex = 0;
  const int myInputIndex = nextInputIndex++;
  int orientIndex = 0;

  // isblinker without an explicit period: find the oscillator's period
  // (stepping under the current rule) so every phase gets tried.
  if (input.isBlinker && input.period == 1) {
    LifeState start = LifeState::Parse(rle, input.centerX, input.centerY);
    LifeState cur = start;
    unsigned found = 0;
    for (unsigned p = 1; p <= 32 && !found; p++) {
      cur.Step();
      if ((cur ^ start).IsEmpty())
        found = p;
    }
    if (found) {
      input.period = found;
      std::cout << "isblinker: detected period " << found << " for " << input.rle << std::endl;
    } else {
      std::cout << "Warning: isblinker catalyst " << input.rle
                << " is not an oscillator of period <= 32 in this rule" << std::endl;
    }
  }

  for (auto &tran : trans) {
    orientIndex++;
    LifeState pat = LifeState::Parse(rle, input.centerX, input.centerY, tran);

    for (unsigned i = 0; i < input.period; i++) {
      CatalystData result;

      result.name = input.rle;
      result.inputIndex = myInputIndex;
      result.orient = orientIndex;
      result.state = pat;

      result.target = LifeTarget(pat);

      result.reactionMask = pat.BigZOI();
      result.reactionMask.Transform(Rotate180OddBoth);
      result.reactionMask.RecalculateMinMax();

      result.maxDisappear = input.maxDisappear;

      for (unsigned k = 0; k < input.forbiddenRLE.size(); k++) {
        result.forbidden.push_back(LifeTarget::Parse(input.forbiddenRLE[k].c_str(),
                                                     input.forbiddenXY[k].first,
                                                     input.forbiddenXY[k].second, tran));
      }

      if (input.requiredRLE != "") {
        result.required = LifeState::Parse(input.requiredRLE.c_str(),
                                           input.requiredXY.first,
                                           input.requiredXY.second, tran);
      }

      if (input.antirequiredRLE != "") {
        result.antirequired = LifeState::Parse(input.antirequiredRLE.c_str(),
                                               input.antirequiredXY.first,
                                               input.antirequiredXY.second, tran);
      }

      if (input.locusRLE != "") {
        result.hasLocus = true;
        result.locus = LifeState::Parse(input.locusRLE.c_str(),
                                        input.locusXY.first,
                                        input.locusXY.second, tran);
      } else {
        result.hasLocus = false;

        // Use the envelope
        LifeState tmp = pat;
        for (unsigned j = 0; j < input.period; j++) {
          result.locus |= tmp;
          tmp.Step();
        }
      }

      result.locusReactionMask = result.locus.BigZOI();
      result.locusReactionMask.Transform(Rotate180OddBoth);
      result.locusReactionMask.RecalculateMinMax();

      // This is a little janky, it should probably be per phase
      result.locusAvoidMask = result.reactionMask;
      result.locusAvoidMask &= ~result.locusReactionMask;
      result.locusAvoidMask.RecalculateMinMax();

      LifeState tmp = pat;
      for (unsigned j = 0; j < input.period; j++) {
        result.phases.push_back(tmp);

        LifeState phasemask = tmp.BigZOI();
        phasemask.Transform(Rotate180OddBoth);
        phasemask &= result.locusReactionMask;
        phasemask.RecalculateMinMax();
        result.phaseReactionMask.push_back(phasemask);

        tmp.Step();
      }

      result.transparent = input.transparent;
      result.maxUses = input.maxUses;
      for (size_t z = 0; z < input.zoneRLE.size(); z++)
        result.zone |= LifeState::Parse(input.zoneRLE[z].c_str(), input.zoneXY[z].first, input.zoneXY[z].second);
      result.period = input.period;
      result.mustInclude = input.mustInclude;
      result.checkRecovery = input.checkRecovery;
      result.sacrificial = input.sacrificial;
      result.isBlinker = input.isBlinker;

      results.push_back(result);

      pat.Step();
    }
  }
  return results;
}

// How a search node's state evolves if no further catalyst is placed. A child
// node differs from its parent only in the columns its new catalyst has
// influenced (at most one more column on each side per generation), so it
// steps only those columns and takes the others from the parent's trajectory.
// After each step the set is narrowed to the columns that really differ; once
// it covers most of the pattern, plain stepping is cheaper and is used instead.
struct Trajectory {
  Trajectory *parent = nullptr;
  row_t diff = 0;       // columns where the current state may differ from the parent's
  unsigned start = 0;   // generation of states[0] (nodes that store their states)
  unsigned limit = 0;   // the trajectory is available up to this generation
  std::vector<LifeState> states;
  std::vector<row_t> changed; // changed[k]: columns of states[k] that differ from states[k - 1]

  const LifeState &At(unsigned gen) {
    while (start + states.size() <= gen) {
      LifeState n = states.back();
      Advance(n);
      row_t ch = 0;
      const LifeState &prev = states.back();
      for (int i = 0; i < N; i++)
        ch |= (row_t)(n.state[i] != prev.state[i]) << i;
      states.push_back(n);
      changed.push_back(ch);
    }
    return states[gen - start];
  }

  // `cur` (this node's state at generation t) -> generation t + 1, in place
  void Advance(LifeState &cur) {
    const unsigned t = cur.gen;
    if (parent && t + 1 <= parent->limit) {
      const row_t cone = diff | RotateLeft(diff) | RotateRight(diff);
      const row_t nz = NonEmptyColumns(cur.state);
      const row_t active = nz | RotateLeft(nz) | RotateRight(nz);
      const row_t live = cone & active;
      if (10 * RowPop(live) <= 9 * RowPop(active)) {
        const LifeState &pn = parent->At(t + 1);
        const row_t pch = parent->changed[t + 1 - parent->start];
        row_t buf[N];
        cur.StepColumnsInto(live, buf);
        for (row_t b = pch & ~cone; b; b &= b - 1) {
          int i = RowCtz(b);
          cur.state[i] = pn.state[i];
        }
        row_t nd = 0;
        for (row_t b = cone; b; b &= b - 1) {
          int i = RowCtz(b);
          row_t v = ((live >> i) & 1) ? buf[i] : 0;
          cur.state[i] = v;
          nd |= (row_t)(v != pn.state[i]) << i;
        }
        diff = nd;
        cur.gen = t + 1;
        cur.min = 0;
        cur.max = N - 1;
        return;
      }
      parent = nullptr; // from here on the difference is too large to be worth it
    }
    cur.Step();
    cur.min = 0;
    cur.max = N - 1;
  }
};

struct Configuration {
  unsigned count;
  unsigned transparentCount;
  unsigned mustIncludeCount;
  std::array<int, MAX_CATALYSTS> curx;
  std::array<int, MAX_CATALYSTS> cury;
  std::array<int, MAX_CATALYSTS> curs;
  LifeState state;
  LifeState catalystsState;
  LifeState startingCatalysts;
  // Bait tracking (inherited by child configurations)
  std::array<unsigned, MAX_BAITS> baitAbsent{};  // generations the bait was missing
  std::array<unsigned, MAX_BAITS> baitStreak{};  // consecutive generations it has been present
  // True if catalystsState is a fixed point (all catalysts are still lifes that do
  // not interact), so stepping it would be a no-op. Recomputed at every placement.
  bool catalystsStill{true};
};

// Fix a, what positions of b causes a collision?
LifeState CollisionMask(const LifeState &a, const LifeState &b, unsigned period) {
  unsigned popsum = a.GetPop() + b.GetPop();

  LifeState mask;
  for (unsigned x = 0; x < N; x++) {
    for (unsigned y = 0; y < N; y++) {
      LifeState state = a;
      state.Join(b, x, y);

      // No overlaps allowed
      if (!(state.GetPop() == popsum)) {
        mask.Set(x, y);
        continue;
      }

      for(unsigned i = 0; i < period; i++)
        state.Step();

      if (!state.Contains(a) || !state.Contains(b, x, y) || state.GetPop() != popsum) {
        mask.Set(x, y);
      }

    }
  }

  return mask;
}

unsigned gcd(unsigned a, unsigned b) {
    while(a != b)
    {
        if(a > b)
            a -= b;
        else
            b -= a;
    }
    return a;
}

// FNV-1a, used to key the on-disk mask caches.
static uint64_t Fnv1a(const std::string &str) {
  uint64_t h = 1469598103934665603ULL;
  for (unsigned char c : str) {
    h ^= c;
    h *= 1099511628211ULL;
  }
  return h;
}

// Exact identity of a LifeState's cells (GetHash() is a cheap checksum that
// different shapes can share, so it must not be used as a cache key).
static uint64_t StateKey(const LifeState &st) {
  uint64_t h = 1469598103934665603ULL;
  for (int i = 0; i < N; i++) {
    row_t w = st.state[i];
    for (int b = 0; b < (int)sizeof(row_t); b++) {
      h ^= (w >> (8 * b)) & 0xFF;
      h *= 1099511628211ULL;
    }
  }
  return h;
}

// Writes a cache file via a temporary name + rename, so an interrupted run never
// leaves a truncated file behind that a later run would read as valid masks.
static void WriteCacheAtomically(const std::string &fname, const char *data, size_t bytes) {
  std::string tmp = fname + ".tmp" + std::to_string((unsigned long long)std::hash<std::thread::id>()(std::this_thread::get_id()) % 100000);
  {
    std::ofstream out(tmp.c_str(), std::ofstream::binary);
    if (!out.good())
      return; // the cache is optional
    out.write(data, (std::streamsize)bytes);
    if (!out.good()) {
      out.close();
      std::remove(tmp.c_str());
      return;
    }
  }
  if (std::rename(tmp.c_str(), fname.c_str()) != 0)
    std::remove(tmp.c_str());
}

// The collision masks depend on the rule and on the lcm of the two periods, not
// just on the two shapes, so both are part of the cache key. (Old caches named
// "maskraw-<hash>-<hash>" / "maskpack-<hash>" are ignored.)
LifeState LoadCollisionMask(const CatalystData &a, const CatalystData &b) {
  unsigned lcmPeriod = (a.period * b.period) / gcd(a.period, b.period);
  std::stringstream ss;
  ss << "masks/maskraw3-" << (N == 64 ? "" : "t" + std::to_string(N) + "-") << std::hex << Fnv1a(RuleTableKey()) << "-" << StateKey(a.state) << "-"
     << StateKey(b.state) << std::dec << "-p" << lcmPeriod;
  std::string fname = ss.str();
  // The per-pair cache is opt-in (as before): it is only used if a directory
  // named "masks" already exists, because big catalyst lists create one small
  // file per pair. The whole-list "maskpack2-*" file is always written.

  std::ifstream infile;
  infile.open(fname.c_str(), std::ios::binary);
  if (!infile.good()) {
    LifeState mask = CollisionMask(a.state, b.state, lcmPeriod);
    WriteCacheAtomically(fname, (const char *)mask.state, N * sizeof(row_t)); // silently skipped without masks/
    return mask;
  } else {
    LifeState result;
    infile.read((char*)result.state, N * sizeof(row_t));
    if (!infile) { // truncated/corrupt cache file: recompute
      infile.close();
      LifeState mask = CollisionMask(a.state, b.state, lcmPeriod);
      WriteCacheAtomically(fname, (const char *)mask.state, N * sizeof(row_t));
      return mask;
    }
    return result;
  }
}

// RLE body of a grid indexed [x][y]. Cell values are states: 0 = dead, 1 = alive
// (written b/o), or with `multiState` any state 0..24 (written . A B C ...),
// as used by multi-state rules such as LifeHistory.
std::string GetRLE(const std::vector<std::vector<uint8_t>> &life2d, bool multiState = false) {
  if (life2d.empty())
    return "";

  if (life2d[0].empty())
    return "";

  auto symbol = [multiState](uint8_t v) -> char {
    if (multiState)
      return v == 0 ? '.' : (char)('A' + v - 1);
    return v ? 'o' : 'b';
  };

  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < life2d[0].size(); j++) {
    uint8_t last_val = life2d[0][j];
    unsigned run_count = 0;

    for (const auto & i : life2d) {
      uint8_t val = i[j];

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val != last_val) {
        if (run_count > 1)
          result << run_count;

        result << symbol(last_val);

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val) {
      if (run_count > 1)
        result << run_count;

      result << symbol(last_val);

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}

// What the output needs to know about a solution besides its cells.
struct SolutionInfo {
  LifeState cats;          // the catalysts at generation 0 (with symmetric images)
  std::string desc;        // the output-list line
  std::vector<int> usage;  // cat lines used (0-based input indices)
  unsigned catCount = 0;   // catalysts in it (fixed ones included)
  bool edge = false;       // its reaction reaches the torus edge
};

class SearchResult {
public:
  SolutionInfo info;
  unsigned junk = 0;       // non-catalyst cells of the outcome at the end
  // Saved for the report
  LifeState init;
  // Every cell the result touches while it runs (init through the container's
  // horizon). Frames in the output are sized by this, so results drawn side by
  // side do not run into each other when the output is simulated in Golly.
  LifeState extent;

  // iters state in form of integers
  // std::vector<int> params;
  unsigned maxGenSurvive;
  unsigned firstGenSurvive;
  // Cells where the goal (filter / match pattern) was found, at the generation it was found.
  LifeState goal;

  SearchResult(LifeState &initState, const Configuration &conf,
               unsigned firstGenSurviveIn, unsigned genSurvive) {
    init = initState;

    maxGenSurvive = genSurvive;
    firstGenSurvive = firstGenSurviveIn;
  }

  void Print() {
    std::cout << "start:" << firstGenSurvive;
    std::cout << ", finish:" << maxGenSurvive << ", params: ";

    // for (int param : params)
    //   std::cout << param << ",";

    std::cout << std::endl;
  }
};

class Category {
private:
  unsigned catDelta;
  int maxgen;
  uint64_t hash;

public:
  LifeState categoryKey;
  std::vector<SearchResult> results;

  Category(LifeState &catalystRemoved, SearchResult &firstResult,
           unsigned catDeltaIn, unsigned maxGen) {
    categoryKey = catalystRemoved;
    results.push_back(firstResult);
    catDelta = catDeltaIn;
    maxgen = maxGen;

    LifeState temp = categoryKey;
    temp.Step(maxgen - temp.gen);
    hash = temp.GetHash();
  }

  void Add(SearchResult &result) { results.push_back(result); }

  bool BelongsTo(LifeState &test, const uint64_t &testHash) {
    if (testHash != hash)
      return false;

    LifeState tempCat = categoryKey;
    LifeState tempTest = test;

    if (tempTest.gen > tempCat.gen)
      tempCat.Step(tempTest.gen - tempCat.gen);
    else if (tempTest.gen < tempCat.gen)
      tempTest.Step(tempCat.gen - tempTest.gen);

    for (unsigned i = 0; i < catDelta; i++) {
      if (tempTest == tempCat)
        return true;

      tempCat.Step();
      tempTest.Step();
    }

    return false;
  }

  static bool CompareSearchResult(SearchResult &a, SearchResult &b) {
    return (a.maxGenSurvive - a.firstGenSurvive) >
           (b.maxGenSurvive - b.firstGenSurvive);
  }

  void Sort() {
    std::sort(results.begin(), results.end(), Category::CompareSearchResult);
  }

  void Print() {
    for (auto & result : results)
      result.Print();
  }

  // Lays the results out in one row. Every result is drawn in a common frame:
  // the bounding box of everything the results in this category touch while
  // they run (not just their starting cells), plus `gap` empty cells. That keeps
  // the results aligned and lets the output be run in Golly without
  // neighbouring results colliding (escaping spaceships excepted).
  // With `historyPat` the row is drawn in LifeHistory states: catalysts state 1
  // (green), the starting pattern `historyPat` state 3 (white, leaves red
  // state-4 cells where it dies), and the cells where the goal was found
  // state 4 (red) if they start dead.
  std::string RLE(int maxCatSize, int gap, const LifeState *historyPat = nullptr, unsigned *rowHeight = nullptr) {
    unsigned howmany = results.size();

    if (maxCatSize != -1)
      howmany = std::min(howmany, (unsigned)maxCatSize);

    LifeState envelope;
    for (unsigned l = 0; l < howmany; l++)
      envelope |= historyPat ? results[l].extent | results[l].goal : results[l].extent;

    int minI = N, maxI = -1, minJ = N, maxJ = -1;
    for (int j = 0; j < N; j++)
      for (int i = 0; i < N; i++)
        if (envelope.GetCell(i - N / 2, j - N / 2) == 1) {
          minI = std::min(minI, i);
          maxI = std::max(maxI, i);
          minJ = std::min(minJ, j);
          maxJ = std::max(maxJ, j);
        }

    if (maxI < 0) { // nothing alive; keep a valid (empty) frame
      minI = minJ = 0;
      maxI = maxJ = 0;
    }

    const unsigned frameW = maxI - minI + 1;
    const unsigned frameH = maxJ - minJ + 1;
    const unsigned Dist = frameW + gap;

    unsigned width = Dist * howmany;
    unsigned height = frameH + gap;
    if (rowHeight)
      *rowHeight = height;

    std::vector<std::vector<uint8_t>> vec(width, std::vector<uint8_t>(height));

    for (unsigned l = 0; l < howmany; l++)
      for (int j = minJ; j <= maxJ; j++)
        for (int i = minI; i <= maxI; i++) {
          const int x = i - N / 2, y = j - N / 2;
          uint8_t v = results[l].init.GetCell(x, y);
          if (historyPat) {
            if (v)
              v = historyPat->GetCell(x, y) ? 3 : 1;
            else if (results[l].goal.GetCell(x, y))
              v = 4;
          }
          vec[Dist * l + (i - minI)][j - minJ] = v;
        }

    return GetRLE(vec, historyPat != nullptr);
  }
};

class CategoryContainer {
public:
  std::vector<Category *> categories;
  unsigned catDelta;
  unsigned maxgen;

  // Duplicate removal. A solution is a duplicate if its starting state is
  // identical to one already stored, or is the image of one under a symmetry
  // that leaves the active pattern (and also-required cells) unchanged.
  bool dedupEnabled = true;
  // Symmetries of the problem: transform, then shift by (dx, dy).
  struct DedupSym { SymmetryTransform t; int dx, dy; };
  std::vector<DedupSym> dedupTransforms{{Identity, 0, 0}};
  std::unordered_map<uint64_t, std::vector<LifeState>> seenInits;
  unsigned duplicatesDropped = 0;

  explicit CategoryContainer(unsigned maxGen) {
    catDelta = 14;
    maxgen = maxGen + catDelta;
  }

  CategoryContainer(unsigned cats, unsigned maxGen) {
    catDelta = cats;
    maxgen = maxGen + catDelta;
  }

  // Every cell `init` touches in its first `gens` generations
  static LifeState Extent(const LifeState &init, unsigned gens) {
    LifeState extent = init, cur = init;
    for (unsigned g = 0; g < gens; g++) {
      cur.Step();
      extent |= cur;
    }
    return extent;
  }

  bool IsDuplicate(const LifeState &init) const {
    if (!dedupEnabled)
      return false;
    for (const DedupSym &sym : dedupTransforms) {
      LifeState image = init;
      if (sym.t != Identity)
        image.Transform(sym.t);
      if (sym.dx != 0 || sym.dy != 0)
        image.Move(sym.dx, sym.dy);
      auto it = seenInits.find(image.GetHash());
      if (it == seenInits.end())
        continue;
      for (const LifeState &seen : it->second)
        if (seen == image)
          return true;
    }
    return false;
  }

  // Returns false (and stores nothing) if the solution is a duplicate.
  bool Add(LifeState &init, const LifeState &afterCatalyst,
           const LifeState &catalysts, const Configuration &conf,
           unsigned firstGenSurvive,
           unsigned genSurvive, const LifeState &goal = LifeState(), const LifeState *knownExtent = nullptr,
           const SolutionInfo &info = SolutionInfo()) {
    if (dedupEnabled) {
      if (IsDuplicate(init)) {
        duplicatesDropped++;
        return false;
      }
      seenInits[init.GetHash()].push_back(init);
    }

    LifeState extent = knownExtent ? *knownExtent : Extent(init, maxgen);

    LifeState result = afterCatalyst ^ catalysts;
    result.gen = firstGenSurvive;

    result.Step(maxgen - result.gen);
    uint64_t hash = result.GetHash();
    const unsigned junk = result.GetPop();

    result = afterCatalyst ^ catalysts;
    result.gen = firstGenSurvive;

    for (auto & category: categories) {
      if (category->BelongsTo(result, hash)) {
          SearchResult r(init, conf, firstGenSurvive, genSurvive);
          r.extent = extent;
          r.goal = goal;
          r.info = info;
          r.junk = junk;
          category->Add(r);
          return true;
      }
    }

    LifeState categoryKey = afterCatalyst ^ catalysts;
    categoryKey.gen = firstGenSurvive;

    SearchResult r(init, conf, firstGenSurvive, genSurvive);
    r.extent = extent;
    r.goal = goal;
    r.info = info;
    r.junk = junk;
    categories.push_back(new Category(categoryKey, r, catDelta, maxgen));
    return true;
  }

  void Sort() {
    for (auto & category: categories)
      category->Sort();
  }

  // Orders the results inside each row and then the rows, by the comma-separated
  // keys of sort-by: found (as found), recovery (earliest first), cells (fewest
  // catalyst cells), catalysts (fewest catalysts), junk (smallest outcome),
  // frequency (rows with the most results first).
  void SortResults(const std::string &spec) {
    std::vector<std::string> keys;
    split(spec, ',', keys);
    if (keys.empty() || (keys.size() == 1 && keys[0] == "found"))
      return;
    auto cmpRes = [&](const SearchResult &a, const SearchResult &b) {
      for (const std::string &k : keys) {
        long d = 0;
        if (k == "recovery") d = (long)a.firstGenSurvive - (long)b.firstGenSurvive;
        else if (k == "cells") d = (long)a.info.cats.GetPop() - (long)b.info.cats.GetPop();
        else if (k == "catalysts") d = (long)a.info.catCount - (long)b.info.catCount;
        else if (k == "junk") d = (long)a.junk - (long)b.junk;
        if (d) return d < 0;
      }
      return false;
    };
    for (auto &c : categories)
      std::stable_sort(c->results.begin(), c->results.end(), cmpRes);
    std::stable_sort(categories.begin(), categories.end(), [&](const Category *a, const Category *b) {
      for (const std::string &k : keys) {
        if (k == "frequency" && a->results.size() != b->results.size())
          return a->results.size() > b->results.size();
        if (k != "frequency" && k != "found") {
          if (cmpRes(a->results[0], b->results[0])) return true;
          if (cmpRes(b->results[0], a->results[0])) return false;
        }
      }
      return false;
    });
  }

  // Machine-readable list of every result, read back by --merge.
  std::string DataFile(const std::string &header) const {
    std::stringstream ss;
    ss << "# CatForce results data 1\n" << header;
    ss << "# first-gen<TAB>catalysts rle x y<TAB>goal rle x y<TAB>cat lines<TAB>catalysts<TAB>edge<TAB>description\n";
    for (const Category *c : categories)
      for (const SearchResult &r : c->results) {
        int x0, y0, w, h, gx, gy;
        std::string cats = BoundingBoxRle(r.info.cats, x0, y0, w, h);
        std::string goal = r.goal.IsEmpty() ? "-" : BoundingBoxRle(r.goal, gx, gy, w, h);
        ss << r.firstGenSurvive << "\t" << cats << " " << x0 << " " << y0 << "\t";
        if (r.goal.IsEmpty()) ss << "- 0 0"; else ss << goal << " " << gx << " " << gy;
        ss << "\t";
        if (r.info.usage.empty()) ss << "-";
        for (size_t u = 0; u < r.info.usage.size(); u++) ss << (u ? "," : "") << r.info.usage[u];
        ss << "\t" << r.info.catCount << "\t" << (r.info.edge ? 1 : 0) << "\t" << r.info.desc << "\n";
      }
    return ss.str();
  }

  void Print() {
    for (auto & category: categories)
      category->Print();
  }

  // With `summary`, also one "#C row ..." line per row (where it starts and what is in it).
  std::string CategoriesRLE(int maxCatSize, int gap, const LifeState *historyPat = nullptr,
                            std::string *summary = nullptr) {
    std::stringstream ss, sum;
    unsigned y = 0, k = 0;
    for (auto & category: categories) {
      unsigned h = 0;
      ss << category->RLE(maxCatSize, gap, historyPat, &h);
      if (summary && k < 500) {
        const SearchResult &r = category->results[0];
        sum << "#C row " << (k + 1) << " at y = " << y << ": " << category->results.size() << " result"
            << (category->results.size() > 1 ? "s" : "") << "; first: recovered from gen " << r.firstGenSurvive << ", "
            << r.info.catCount << " catalyst" << (r.info.catCount == 1 ? "" : "s") << " (" << r.info.cats.GetPop()
            << " cells), outcome " << r.junk << " cells";
        if (!r.info.usage.empty()) {
          sum << ", cat lines";
          for (size_t u = 0; u < r.info.usage.size(); u++)
            sum << (u ? "," : " ") << (r.info.usage[u] + 1);
        }
        sum << (r.info.edge ? ", reaches the torus edge" : "") << "\n";
      } else if (summary && k == 500) {
        sum << "#C (more rows not listed)\n";
      }
      y += h;
      k++;
    }
    if (summary)
      *summary = sum.str();
    ss << "!";

    return ss.str();
  }

  // One square cell per distinct outcome (first solution of each row), cells laid out in a square grid.
  std::string UniqueRLE(int gap) {
    // Each cell shows the starting pattern (pattern + catalysts) of the first solution of a row,
    // cropped to its own bounding box, so the reaction envelope does not blow the cells up.
    unsigned n = categories.size(), cols = 1;
    while (cols * cols < n) cols++;
    std::vector<std::array<int, 4>> box(n);
    int side = 1;
    for (unsigned k = 0; k < n; k++) {
      int minI = N, maxI = -1, minJ = N, maxJ = -1;
      for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
          if (categories[k]->results[0].init.GetCell(i - N / 2, j - N / 2)) {
            minI = std::min(minI, i); maxI = std::max(maxI, i);
            minJ = std::min(minJ, j); maxJ = std::max(maxJ, j);
          }
      if (maxI < 0) minI = minJ = maxI = maxJ = 0;
      box[k] = {minI, maxI, minJ, maxJ};
      side = std::max(side, std::max(maxI - minI, maxJ - minJ) + 1);
    }
    side += gap;
    unsigned rows = n ? (n + cols - 1) / cols : 1;
    std::vector<std::vector<uint8_t>> vec(side * cols, std::vector<uint8_t>(side * rows));
    for (unsigned k = 0; k < n; k++)
      for (int j = box[k][2]; j <= box[k][3]; j++)
        for (int i = box[k][0]; i <= box[k][1]; i++)
          vec[side * (k % cols) + (i - box[k][0])][side * (k / cols) + (j - box[k][2])] =
              categories[k]->results[0].init.GetCell(i - N / 2, j - N / 2);
    return GetRLE(vec) + "!";
  }
};


class CatalystSearcher {
public:
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  mutable std::mutex mtx;                       // guards results, counters, reporting
  std::atomic<bool> stopFlag{false};
  std::string stopReason;
  mutable std::chrono::steady_clock::time_point lastReportTime{};
  std::vector<std::string> solutionLines;       // for output-list
  bool searchFinished = false;
  std::vector<std::string> inputNames;          // cat line RLEs, by input index
  std::vector<unsigned> catalystUsage;          // solutions using each cat line

  double ElapsedSeconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  }

  // Polled by the search; true once Ctrl-C, max-time or max-results hit.
  bool ShouldStop() {
    if (gInterrupted && !stopFlag.exchange(true)) {
      std::lock_guard<std::mutex> lk(mtx);
      stopReason = "interrupted (Ctrl-C)";
      std::cout << std::endl << "Stopping and saving... (press Ctrl-C again to quit at once)" << std::endl;
    }
    if (stopFlag.load(std::memory_order_relaxed))
      return true;
    if (params.maxTimeSec > 0) {
      static thread_local unsigned tick = 0;
      if ((++tick & 1023) == 0 && ElapsedSeconds() > params.maxTimeSec) {
        if (!stopFlag.exchange(true)) {
          std::lock_guard<std::mutex> lk(mtx);
          stopReason = "max-time reached";
        }
        return true;
      }
    }
    return false;
  }
  SearchParams params;
  LifeState pat;
  LifeState alsoRequired;
  std::vector<TargetSet> baitTargets;      // per bait: accepted phases of the (symmetry-joined) bait
  std::vector<unsigned> baitSurvive;
  std::vector<unsigned> baitAbsentMin;
  std::vector<CatalystData> catalysts;
  std::vector<LifeTarget> targetFilter;
  std::vector<LifeState> catalystCollisionMasks;
  // Latest generation at which placing each catalyst can still lead to a
  // solution (see ComputeLatestPlacementGens), and the latest over all of them.
  std::vector<int> latestPlacementGen;
  int latestPlacementGenAny = 0;
  // fixed-cat entries: their index in `catalysts` and their position
  std::vector<unsigned> fixedIdx;
  std::vector<std::pair<int, int>> fixedPos;
  unsigned rootCount = 0;  // catalysts in place before the search starts (the fixed ones)
  unsigned totalCats = 0;  // rootCount + num-catalyst
  LifeState catalystZone, noCatalystZone, keepEmpty;
  std::atomic<unsigned> edgeFlagged{0}, edgeDropped{0};

  long long idx{};
  unsigned found{};
  unsigned fullfound{};
  long long total{};
  unsigned short int counter{};

  CategoryContainer *categoryContainer{};
  CategoryContainer *fullCategoryContainer{};

  bool hasFilter{};
  bool hasMustInclude{};
  bool reportAll{};

  unsigned filterMaxGen{};

  // Key for the on-disk mask pack. Must change whenever any mask could change:
  // the rule, the catalysts *in order* (the pack is indexed by position), their
  // periods, and the options that make LoadMasks skip some pairs.
  uint64_t AllCatalystsHash() const {
    std::stringstream ss;
    // (only the exact pair-skipping conditions used in LoadMasks matter, not num-catalyst itself)
    ss << RuleTableKey() << (N == 64 ? "" : "|torus" + std::to_string(N)) << "|fixed" << fixedIdx.size() << "|skipMust" << (params.numCatalysts == 2 && hasMustInclude)
       << "|skipTransp" << (params.numTransparent == 1);
    for (auto &cat : catalysts)
      ss << "|" << StateKey(cat.state) << "," << cat.period << "," << cat.mustInclude << "," << cat.transparent << "," << cat.fixed;
    return Fnv1a(ss.str());
  }

  void LoadMasks() {
    if (totalCats == 1)
      return;

    catalystCollisionMasks = std::vector<LifeState>(catalysts.size() * catalysts.size());

    std::stringstream ss;
    ss << "maskpack3-" << std::hex << AllCatalystsHash();
    std::string fname = ss.str();

    // See if the pack exists
    std::ifstream infile;
    infile.open(fname.c_str(), std::ios::binary);
    if (infile.good()) {
      for (unsigned s = 0; s < catalysts.size(); s++) {
        for (unsigned t = 0; t < catalysts.size(); t++) {
          infile.read((char*)catalystCollisionMasks[s * catalysts.size() + t].state, N * sizeof(row_t));
          catalystCollisionMasks[s * catalysts.size() + t].RecalculateMinMax();
        }
      }
      // A short or damaged file must not be trusted: recompute instead.
      if (infile.good() && infile.peek() == std::ifstream::traits_type::eof())
        return;
      std::cout << "Mask cache " << fname << " is damaged; recomputing it." << std::endl;
      catalystCollisionMasks.assign(catalysts.size() * catalysts.size(), LifeState());
    }
    infile.close();

    // If not, load or generate the masks
    for (unsigned s = 0; s < catalysts.size(); s++) {
      for (unsigned t = 0; t < catalysts.size(); t++) {
        if (catalysts[t].fixed)
          continue; // never placed next to anything
        if (!catalysts[s].fixed) {
          if (params.numCatalysts == 2 && hasMustInclude &&
              !catalysts[s].mustInclude && !catalysts[t].mustInclude)
            continue;

          if (params.numTransparent == 1 && catalysts[s].transparent &&
              catalysts[t].transparent)
            continue;
        }

        catalystCollisionMasks[s * catalysts.size() + t] = LoadCollisionMask(catalysts[s], catalysts[t]);
        catalystCollisionMasks[s * catalysts.size() + t].RecalculateMinMax();
      }
    }

    // Save the pack for next time (atomically, see WriteCacheAtomically)
    std::vector<row_t> buffer;
    buffer.reserve(catalysts.size() * catalysts.size() * N);
    for (unsigned s = 0; s < catalysts.size(); s++)
      for (unsigned t = 0; t < catalysts.size(); t++)
        buffer.insert(buffer.end(), catalystCollisionMasks[s * catalysts.size() + t].state,
                      catalystCollisionMasks[s * catalysts.size() + t].state + N);
    WriteCacheAtomically(fname, (const char *)buffer.data(), buffer.size() * sizeof(row_t));
  }

  void Init(const char *inputFile) {
    begin = std::chrono::steady_clock::now();

    std::vector<CatalystInput> inputcats;
    ReadParams(inputFile, inputcats, params);
    params.inputFile = inputFile;

    for (auto &input : inputcats) {
      inputNames.push_back(input.rle);
      std::vector<CatalystData> newcats = CatalystData::FromInput(input);
      catalysts.insert(catalysts.end(), newcats.begin(), newcats.end());
    }
    // fixed-cat: one entry each (phase 0, as drawn), placed by Search()
    for (const auto &f : params.fixedCats) {
      std::vector<std::string> rest = splitwhitespace(f.rest);
      std::string line = "cat " + f.rle + " " + rest[0] + " 0 0 .";
      for (size_t k = 1; k < rest.size(); k++)
        line += " " + rest[k];
      CatalystInput in(line);
      CatalystData d = CatalystData::FromInput(in)[0];
      d.fixed = true;
      d.inputIndex = -1;
      fixedIdx.push_back(catalysts.size());
      fixedPos.push_back({f.x, f.y});
      catalysts.push_back(d);
    }
    rootCount = fixedIdx.size();
    totalCats = params.numCatalysts + rootCount;
    minTotal = params.minCatalysts + rootCount;
    if (totalCats > (unsigned)MAX_CATALYSTS)
      Die("num-catalyst plus fixed-cat lines may be at most " + std::to_string(MAX_CATALYSTS));

    auto unionOf = [](const std::vector<SearchParams::PlacedRle> &v) {
      LifeState r;
      for (const auto &p : v)
        r |= LifeState::Parse(p.rle.c_str(), p.x, p.y);
      return r;
    };
    catalystZone = unionOf(params.catalystZones);
    noCatalystZone = unionOf(params.noCatalystZones);
    keepEmpty = unionOf(params.keepEmpty);

    hasMustInclude = false;
    for (auto &cat : catalysts) {
      hasMustInclude = hasMustInclude || cat.mustInclude;
    }

    if (!hasMustInclude) {
      for (auto &cat : catalysts) {
        cat.mustInclude = true;
      }
    }

    catalystUsage.assign(inputNames.size(), 0);
    ComputeLatestPlacementGens();
    pat = LifeState::Parse(params.pat.c_str(), params.xPat, params.yPat);
    SetupBaits();
    categoryContainer = new CategoryContainer(params.maxGen);
    fullCategoryContainer = new CategoryContainer(params.maxGen);
    // (also-required must be parsed before ConfigureDuplicateRemoval, which
    // only accepts symmetries that also map the also-required cells onto themselves)
    alsoRequired = LifeState::Parse(params.alsoRequired.c_str(), params.alsoRequiredXY.first, params.alsoRequiredXY.second);
    pat |= alsoRequired; // also-required cells are part of the starting pattern
    ConfigureDuplicateRemoval();

    for (unsigned i = 0; i < params.targetFilter.size(); i++)
      targetFilter.push_back(LifeTarget::Parse(params.targetFilter[i].c_str(),
                                               params.filterdx[i], params.filterdy[i]));

    if (!skipMasks)
      LoadMasks();

    found = 0;

    hasFilter = !params.targetFilter.empty();
    reportAll = params.fullReportFile.length() != 0;

    filterMaxGen = FilterMaxGen();
  }

  // Baits are added to the starting pattern; each must later be present
  // (cells alive, surrounding ring dead) for `survive` consecutive generations
  // after having been missing for at least `absent` generations. Oscillating
  // baits are accepted in any of their phases.
  void SetupBaits() {
    for (const auto &b : params.baits) {
      LifeState bait = LifeState::Parse(b.rle.c_str(), b.x, b.y);
      if (bait.IsEmpty()) {
        std::cout << "Bait " << b.rle << " is empty" << std::endl;
        exit(1);
      }
      pat |= bait;
      pat.RecalculateMinMax(); // |= leaves the cached column range stale

      LifeState joined;
      joined.JoinWSymChain(bait, params.symmetryChain);

      TargetSet ts;
      LifeState cur = joined;
      int period = 0;
      for (int p = 1; p <= 32 && !period; p++) {
        cur.Step();
        if ((cur ^ joined).IsEmpty())
          period = p;
      }
      if (!period) {
        std::cout << "Warning: bait " << b.rle
                  << " is not a still life/oscillator (period <= 32) in this rule, "
                     "so it cannot survive; no results are expected" << std::endl;
        period = 1;
      }
      cur = joined;
      for (int p = 0; p < period; p++) {
        LifeTarget t(cur);
        ts.Add(t);
        cur.Step();
      }
      baitTargets.push_back(ts);
      baitSurvive.push_back(b.survive >= 0 ? (unsigned)b.survive : params.stableInterval);
      baitAbsentMin.push_back((unsigned)std::max(0, b.absent));
      std::cout << "Bait " << b.rle << ": period " << period << ", must be missing >= "
                << baitAbsentMin.back() << " gen(s), then survive " << baitSurvive.back()
                << " gen(s)" << std::endl;
    }
  }

  unsigned FilterMaxGen() {
    unsigned maxGen = params.maxGen;

    for (unsigned j = 0; j < targetFilter.size(); j++) {
      if (params.filterGen[j] >= 0 && params.filterGen[j] > maxGen)
        maxGen = params.filterGen[j];

      if (params.filterGenRange[j].second >= 0 && params.filterGenRange[j].second > maxGen)
        maxGen = params.filterGenRange[j].second;
    }

    return maxGen;
  }

  void Report(const std::string& suffix) {
    std::string temp = params.outputFile;
    params.outputFile = params.outputFile + suffix + std::string(".rle");
    Report();
    params.outputFile = temp;
  }

  mutable std::string lastSavedName;
  mutable unsigned lastSavedFound = ~0u;
  mutable unsigned lastSavedFullFound = ~0u;

  void Report(bool saveFile = true) const {
    unsigned sec = (unsigned)ElapsedSeconds() + 1;
    // unsigned checkPerSecond = idx / (sec * 1000);

    std::cout << "results: " << categoryContainer->categories.size() << "/"
              << found;
    if (params.fullReportFile.length() != 0) {
      std::cout << ", unfiltered: " << fullCategoryContainer->categories.size() << "/"
                << fullfound;
    }
    std::cout << ", now: ";
    PrintTime(sec);
    std::cout << std::setfill(' ');
    if (!searchFinished) {
      const unsigned long long done = rootPlacementsDone.load();
      std::cout << "  (first catalyst at gen " << progressGen.load() << " of " << params.startGen << ".."
                << params.lastGen << ", " << done;
      if (totalRoots > 0) {
        double frac = std::min(1.0, done / (double)totalRoots);
        std::cout << " of " << totalRoots << " first placements, " << std::fixed << std::setprecision(1) << 100 * frac << "%";
        double el = ElapsedSeconds() - countSeconds;
        if (frac >= 0.02 && el >= 10) {
          std::cout << ", roughly ";
          PrintTime((unsigned)(el * (1 - frac) / frac));
          std::cout << std::setfill(' ') << " left";
        }
        std::cout << ")";
      } else {
        std::cout << " first placements done)";
      }
    }
    std::cout << std::endl;
    // std::cout << ", " << std::setprecision(1) << std::fixed << checkPerSecond
    //           << "K/sec" << std::endl;

    // Skip rewriting the (potentially large) output files if nothing was
    // added since the last save to the same file.
    if (saveFile && params.outputFile == lastSavedName && found == lastSavedFound &&
        fullfound == lastSavedFullFound)
      saveFile = false;
    if (saveFile) {
      lastSavedName = params.outputFile;
      lastSavedFound = found;
      lastSavedFullFound = fullfound;
      std::cout << "Saving " << params.outputFile << "... " << std::flush;

      categoryContainer->SortResults(params.sortBy);
      std::string rows;
      std::string body = categoryContainer->CategoriesRLE(params.maxCatSize, params.outputGap, nullptr, &rows);
      SaveTextAtomically(params.outputFile,
                         RleHeaderComments(categoryContainer->categories.size(), found, false) + rows +
                         "x = 0, y = 0, rule = " + RuleForHeader() + "\n" + body);
      SaveTextAtomically(DataFileName(),
                         categoryContainer->DataFile("# input " + params.inputFile + " | rule " + RuleForHeader() +
                                                     " | torus " + std::to_string(N) + " | max-gen " + std::to_string(params.maxGen) +
                                                     " | finished " + (searchFinished ? "yes" : "no") + "\n"));
      std::cout << "Done!" << std::endl;

      if (hasFilter) {
        // Same layout as the output, in LifeHistory states (see Category::RLE)
        LifeState startPattern;
        startPattern.JoinWSymChain(pat, params.symmetryChain);
        SaveTextAtomically(HistoryFileName(),
                           RleHeaderComments(categoryContainer->categories.size(), found, false) +
                           "#C History version: white = starting pattern (turns red where it dies), "
                           "red = where the goal was found, green = catalysts.\n" +
                           "x = 0, y = 0, rule = " + HistoryRuleForHeader() + "\n" +
                           categoryContainer->CategoriesRLE(params.maxCatSize, params.outputGap, &startPattern));
      }

      SaveTextAtomically(UniqueFileName(),
                         "#C CatForce: one square per distinct outcome (" +
                         std::to_string(categoryContainer->categories.size()) + ")\n" +
                         "x = 0, y = 0, rule = " + RuleForHeader() + "\n" +
                         categoryContainer->UniqueRLE(params.outputGap));

      if (!params.outputList.empty()) {
        std::vector<std::string> sorted = solutionLines;
        std::sort(sorted.begin(), sorted.end());
        std::stringstream listFile;
        listFile << "# " << sorted.size() << " solution(s); rule " << RuleForHeader()
                 << "\n# coordinates: torus, origin at the centre, x right, y down\n";
        for (size_t i = 0; i < sorted.size(); i++)
          listFile << (i + 1) << ": " << sorted[i] << "\n";
        SaveTextAtomically(params.outputList, listFile.str());
      }

      if (params.fullReportFile.length() != 0) {
        std::cout << "Saving " << params.fullReportFile << "... " << std::flush;

        SaveTextAtomically(params.fullReportFile,
                           RleHeaderComments(fullCategoryContainer->categories.size(), fullfound, true) +
                           "x = 0, y = 0, rule = " + RuleForHeader() + "\n" +
                           fullCategoryContainer->CategoriesRLE(params.maxCatSize, params.outputGap));
        std::cout << "Done!" << std::endl;
      }
    }
  }

  std::string HistoryFileName() const { return AddSuffix(params.outputFile, "_history"); }
  // results.rle -> results.cfdata (read back by --merge)
  std::string DataFileName() const {
    const std::string &o = params.outputFile;
    size_t slash = o.find_last_of("/\\"), dot = o.find_last_of('.');
    return (dot == std::string::npos || (slash != std::string::npos && dot < slash) ? o : o.substr(0, dot)) + ".cfdata";
  }

  // Golly's [rule]History form of the current rule
  static std::string HistoryRuleForHeader() {
    return (RuleIsIsotropic() ? CurrentRuleCanonical() : currentRuleName) + "History";
  }

  std::string UniqueFileName() const {
    return params.uniqueFile.empty() ? AddSuffix(params.outputFile, "_unique") : params.uniqueFile;
  }

  // "#C" comment lines Golly shows in its info box (File > Show Comments).
  std::string RleHeaderComments(size_t rows, unsigned solutions, bool full) const {
    std::stringstream ss;
    ss << "#C CatForce " << (full ? "full report (filters ignored)" : "results") << ": " << solutions
       << " solution(s) in " << rows << " row(s), one row per distinct outcome.\n";
    ss << "#C input: " << params.inputFile << ", rule " << RuleForHeader() << ", "
       << (params.minCatalysts < params.numCatalysts ? std::to_string(params.minCatalysts) + "-" : std::string())
       << params.numCatalysts << " catalyst(s), max-gen " << params.maxGen;
    if (params.sortBy != "found")
      ss << ", sorted by " << params.sortBy;
    if (params.numParts > 1)
      ss << ", part " << (params.partIndex + 1) << "/" << params.numParts;
    if (N != 64)
      ss << ", torus " << N << "x" << N;
    ss << (searchFinished ? "" : " (search still running when saved)") << "\n";
    if (!params.outputList.empty())
      ss << "#C catalyst positions for every solution: " << params.outputList << "\n";
    return ss.str();
  }

  // Called from the search roots: prints/saves at most about once a second.
  void MaybeReport() {
    if (countOnly)
      return;
    std::lock_guard<std::mutex> lk(mtx);
    auto now = std::chrono::steady_clock::now();
    // every second for the first minute, then every 10 seconds
    auto interval = ElapsedSeconds() < 60 ? std::chrono::seconds(1) : std::chrono::seconds(10);
    if (lastReportTime != std::chrono::steady_clock::time_point{} && now - lastReportTime < interval)
      return;
    lastReportTime = now;
    Report();
  }

  void FinalReport() {
    std::lock_guard<std::mutex> lk(mtx);
    lastSavedFound = ~0u; // always rewrite at the end, so the files say the search finished
    Report();
  }

  static void PrintTime(unsigned sec) {
    unsigned hr = sec / 3600;
    unsigned min = (sec / 60) - hr * 60;
    unsigned secs = sec - 3600 * hr - 60 * min;
    std::cout << std::setfill('0');
    std::cout << hr << ":" << std::setw(2) << min << ":" << std::setw(2) << secs;
 }

  bool HasForbidden(Configuration &conf, unsigned curIter) {
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    for (unsigned i = 0; i <= curIter + 1; i++) {
      for (unsigned j = 0; j < conf.count; j++) {
        for (unsigned k = 0; k < catalysts[conf.curs[j]].forbidden.size(); k++) {
          if (workspace.Contains(catalysts[conf.curs[j]].forbidden[k], conf.curx[j], conf.cury[j]))
            return true;
        }
      }
      workspace.Step();
    }

    return false;
  }

  bool FilterForCurrentGenFail(LifeState &workspace) {
    for (unsigned j = 0; j < targetFilter.size(); j++) {
      if (workspace.gen == params.filterGen[j] &&
          workspace.Contains(targetFilter[j]) == false) {
        return true;
      }
    }

    return false;
  }

  // If `goalCells` is given, it receives the cells where each passed filter's
  // pattern was found (at the generation it was found).
  bool ValidateFilters(Configuration &conf, unsigned successtime, unsigned failuretime,
                       LifeState *goalCells = nullptr) {
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    LifeState catalystState = conf.startingCatalysts;

    unsigned maxMatchingPop;
    if(params.maxJunk != -1)
      maxMatchingPop = workspace.GetPop() + params.maxJunk;
    else
      maxMatchingPop = 10000;

    std::vector<bool> filterPassed(params.filterGen.size(), false);

    unsigned stopTime;
    if (params.stopAfterCatsDestroyed != -1)
      stopTime = failuretime + params.stopAfterCatsDestroyed;
    else
      stopTime = filterMaxGen;

    for (unsigned g = 0; g <= std::min(filterMaxGen, stopTime); g++) {
      for (unsigned k = 0; k < params.filterGen.size(); k++) {
        if (filterPassed[k])
          continue; // No need to check it again.

        bool inSingle = workspace.gen == params.filterGen[k];
        bool inRange = params.filterGen[k] == -1 &&
                       params.filterGenRange[k].first <= workspace.gen &&
                       params.filterGenRange[k].second >= workspace.gen;
        bool shouldCheck = inSingle || (inRange && workspace.gen >= successtime);

        bool succeeded = false;
        LifeState junk;
        LifeState foundCells;

        // See whether there is a match at all
        if (shouldCheck && (params.filterType[k] == ANDFILTER ||
                            params.filterType[k] == ORFILTER)) {
          succeeded = workspace.Contains(targetFilter[k]);
          junk = workspace & ~targetFilter[k].wanted & ~catalystState;
          foundCells = targetFilter[k].wanted;
        }

        if (shouldCheck && (params.filterType[k] == MATCHFILTER)) {
          // if(workspace.GetPop() <= maxMatchingPop) {
            LifeState withoutCatalysts = workspace & ~catalystState;
            for (auto sym : SymmetryGroupFromEnum(StaticSymmetry::D8)) {
              LifeTarget transformed = targetFilter[k];
              transformed.Transform(sym);
              LifeState matches = withoutCatalysts.Match(transformed);
              if(!matches.IsEmpty()) {
                succeeded = true;
                foundCells = matches.Convolve(transformed.wanted) & withoutCatalysts;
                junk = withoutCatalysts & ~foundCells;
                break;
              }
            }
          // }
        }

        if (succeeded && (params.maxJunk == -1 || junk.GetPop() <= params.maxJunk)) {
          filterPassed[k] = true;
          if (goalCells)
            *goalCells |= foundCells;

          // If this was an OR filter, consider all the other OR filters passed
          // too.
          if (params.filterType[k] == ORFILTER) {
            for (unsigned j = 0; j < params.filterGen.size(); j++) {
              if (params.filterType[j] == ORFILTER) {
                filterPassed[j] = true;
              }
            }
          }
        }

        // Bail early
        if (workspace.gen == params.filterGen[k] &&
            params.filterType[k] == ANDFILTER &&
            !workspace.Contains(targetFilter[k])
            )
          return false;
      }

      workspace.Step();
      catalystState.Step();
    }

    for (unsigned k = 0; k < params.filterGen.size(); k++)
      if (!filterPassed[k])
        return false;

    return true;
  }

  // The catalysts of a configuration must not interact with each other. The
  // precomputed collision masks only cover a catalyst against another catalyst
  // as placed; with a global symmetry the mirrored/rotated images are extra
  // catalysts that can overlap or touch the originals. So: evolve all the
  // catalysts alone for a full period and require them to come back unchanged.
  static unsigned CatalystPeriodLcm(const std::vector<CatalystData> &cats, const Configuration &conf) {
    unsigned P = 1;
    for (unsigned j = 0; j < conf.count; j++) {
      unsigned p = cats[conf.curs[j]].period;
      P = P / gcd(P, p) * p;
    }
    return P;
  }

  static bool CatalystsCoexist(const LifeState &catsAtPhase, unsigned P) {
    LifeState cur = catsAtPhase;
    for (unsigned i = 0; i < P; i++)
      cur.Step();
    return cur == catsAtPhase;
  }

  std::atomic<unsigned> interactingRejected{0};
  // Progress of the top-level loop (first catalyst), shown in the status line.
  std::atomic<unsigned> progressGen{0};
  std::atomic<unsigned long long> rootPlacementsDone{0};
  unsigned long long totalRoots = 0;     // first placements of this run (from the counting pass)
  double countSeconds = 0;               // time the counting pass took
  std::atomic<unsigned long long> countedRoots{0};
  bool countOnly = false;                // counting pass: count first placements, do not search them
  unsigned minTotal = 0;                 // rootCount + the fewest catalysts per result
  bool skipMasks = false;                // --merge: no search, so no collision masks

  // A catalyst placed at generation g is first counted as recovered at g (or at
  // the next multiple of its period), and a solution needs it recovered for
  // stable-interval generations in a row by generation max-gen - 1. With period
  // p that takes ceil(stable-interval / p) checks, the last of them at least
  // (ceil(stable-interval / p) - 1) * p generations after g. Placements later
  // than that can never be part of a solution, so they are not tried at all.
  // Sacrificial catalysts need not recover and are never limited.
  void ComputeLatestPlacementGens() {
    latestPlacementGen.assign(catalysts.size(), 0);
    latestPlacementGenAny = -1;
    for (size_t s = 0; s < catalysts.size(); s++) {
      unsigned p = std::max(1u, catalysts[s].period);
      unsigned si = params.stableInterval;
      unsigned need = si == 0 ? 0 : ((si + p - 1) / p - 1) * p;
      latestPlacementGen[s] = catalysts[s].sacrificial ? (int)params.maxGen - 1 : (int)params.maxGen - 1 - (int)need;
      if (!catalysts[s].fixed)
        latestPlacementGenAny = std::max(latestPlacementGenAny, latestPlacementGen[s]);
    }
    if (latestPlacementGenAny < (int)params.startGen)
      std::cout << "Warning: no catalyst can be placed early enough to stay recovered for stable-interval ("
                << params.stableInterval << ") generations before max-gen (" << params.maxGen
                << "); no results are possible. Raise max-gen or lower stable-interval/start-gen." << std::endl;
  }

  // Finds the symmetries that map the whole problem (active pattern, global
  // symmetry images, also-required cells) onto itself. Solutions related by
  // such a symmetry are equivalent when the rule is isotropic.
  // Finds every rotation/reflection (about any centre) that maps the whole
  // problem onto itself: the active pattern with its global-symmetry images and
  // any baits, plus the also-required cells. When the rule is isotropic,
  // solutions related by such a symmetry are the same solution.
  void ConfigureDuplicateRemoval() {
    bool enable = !params.keepDuplicates;
    std::vector<CategoryContainer::DedupSym> syms{{Identity, 0, 0}};

    if (enable && RuleIsIsotropic()) {
      LifeState fullPat;
      fullPat.JoinWSymChain(pat, params.symmetryChain);

      // one live cell of the pattern, used to find the shift that lines a
      // transformed copy up with the original
      int rx = 99, ry = 99;
      for (int x = -N / 2; x < N / 2 && rx == 99; x++)
        for (int y = -N / 2; y < N / 2; y++)
          if (fullPat.GetCell(x, y) == 1) { rx = x; ry = y; break; }

      const SymmetryTransform d4[] = {ReflectAcrossX, ReflectAcrossY, Rotate90, Rotate270,
                                      Rotate180OddBoth, ReflectAcrossYeqX, ReflectAcrossYeqNegX};
      for (SymmetryTransform t : d4) {
        LifeState tp = fullPat;
        tp.Transform(t);
        bool found = false;
        for (int x = -N / 2; x < N / 2 && !found; x++)
          for (int y = -N / 2; y < N / 2 && !found; y++) {
            if (tp.GetCell(x, y) != 1)
              continue;
            int dx = rx - x, dy = ry - y;
            LifeState moved = tp;
            moved.Move(dx, dy);
            if (!(moved == fullPat))
              continue;
            LifeState ar = alsoRequired;
            ar.Transform(t);
            ar.Move(dx, dy);
            if (ar == alsoRequired) {
              syms.push_back({t, dx, dy});
              found = true;
            }
          }
      }
    }

    for (CategoryContainer *c : {categoryContainer, fullCategoryContainer}) {
      c->dedupEnabled = enable;
      c->dedupTransforms = syms;
    }
    if (enable)
      std::cout << "Duplicate removal: on (pattern has " << syms.size()
                << " symmetr" << (syms.size() == 1 ? "y" : "ies")
                << " incl. identity; add \"keep-duplicates yes\" to disable)" << std::endl;
  }


  void ReportSolution(Configuration &conf, unsigned successtime, unsigned failuretime) {
    if (HasForbidden(conf, successtime + 3))
      return;

    if (!CatalystsCoexist(conf.startingCatalysts, CatalystPeriodLcm(catalysts, conf))) {
      interactingRejected++;
      if (!params.keepInteracting)
        return;
    }

    // if reportAll - ignore filters and update fullReport
    if (reportAll) {
      LifeState workspace = conf.startingCatalysts;
      workspace.JoinWSymChain(pat, params.symmetryChain);

      LifeState init = workspace;

      workspace.Step(successtime - params.stableInterval + 2);
      LifeState afterCatalyst = workspace;

      LifeState steppedCatalysts = conf.startingCatalysts;
      steppedCatalysts.Step(successtime - params.stableInterval + 2);

      std::lock_guard<std::mutex> lk(mtx);
      if (fullCategoryContainer->Add(init, afterCatalyst, steppedCatalysts, conf,
                                     successtime - params.stableInterval + 2, 0))
        fullfound++;
    }

    // If has filter validate them;
    LifeState goal;
    if (hasFilter) {
      if (!ValidateFilters(conf, successtime, failuretime, &goal))
        return;
    }

    // If all filters validated update results
    LifeState workspace = conf.startingCatalysts;
    workspace.JoinWSymChain(pat, params.symmetryChain);

    LifeState init = workspace;

    workspace.Step(successtime - params.stableInterval + 2);

    LifeState afterCatalyst = workspace;

    LifeState steppedCatalysts = conf.startingCatalysts;
    steppedCatalysts.Step(successtime - params.stableInterval + 2);

    std::string description = DescribeSolution(conf, successtime);

    // A reaction that reaches the edge of the torus may come back from the
    // other side, which cannot happen on the infinite plane.
    LifeState extent = CategoryContainer::Extent(init, categoryContainer->maxgen);
    bool edge = false;
    if (params.torusEdge != "off") {
      for (int k = 0; k < N && !edge; k++)
        edge = extent.GetCell(-N / 2, k) || extent.GetCell(N / 2 - 1, k) || extent.GetCell(k, -N / 2) ||
               extent.GetCell(k, N / 2 - 1);
      if (edge && params.torusEdge == "drop") {
        edgeDropped++;
        return;
      }
      if (edge)
        description += " | REACHES THE TORUS EDGE (may wrap around)";
    }

    SolutionInfo info;
    info.cats = conf.startingCatalysts;
    info.desc = description;
    info.catCount = conf.count;
    info.edge = edge;
    for (unsigned j = 0; j < conf.count; j++) {
      int idx = catalysts[conf.curs[j]].inputIndex;
      if (idx >= 0 && std::find(info.usage.begin(), info.usage.end(), idx) == info.usage.end())
        info.usage.push_back(idx);
    }

    std::lock_guard<std::mutex> lk(mtx);
    if (AddSolution(init, afterCatalyst, steppedCatalysts, conf, successtime - params.stableInterval + 2, goal, &extent, info)) {

      if (params.maxResults > 0 && found >= params.maxResults && !stopFlag.exchange(true))
        stopReason = "max-results reached";
    }
  }

  // Stores a solution in the filtered results (caller holds mtx). False if it is a duplicate.
  bool AddSolution(LifeState &init, const LifeState &afterCatalyst, const LifeState &steppedCatalysts,
                   const Configuration &conf, unsigned firstGen, const LifeState &goal, const LifeState *extent,
                   const SolutionInfo &info) {
    if (!categoryContainer->Add(init, afterCatalyst, steppedCatalysts, conf, firstGen, 0, goal, extent, info))
      return false;
    if (info.edge)
      edgeFlagged++;
    found++;
    solutionLines.push_back(info.desc);
    for (int idx : info.usage)
      if (idx >= 0 && idx < (int)catalystUsage.size())
        catalystUsage[idx]++;
    return true;
  }

  // One line per solution: where each catalyst sits (oriented RLE + top-left).
  std::string DescribeSolution(const Configuration &conf, unsigned successtime) const {
    std::stringstream ss;
    ss << "done by gen " << successtime;
    if (!params.symmetryChain.empty())
      ss << " (each catalyst also has its symmetric image(s))";
    for (unsigned j = 0; j < conf.count; j++) {
      const CatalystData &c = catalysts[conf.curs[j]];
      LifeState placed = c.state;
      placed.Move(conf.curx[j], conf.cury[j]);
      int x0, y0, w, h;
      std::string rle = BoundingBoxRle(placed, x0, y0, w, h);
      if (c.fixed)
        ss << " | fixed-cat [" << c.name << "] top-left (" << x0 << "," << y0 << ")";
      else
        ss << " | cat line " << (c.inputIndex + 1) << " [" << c.name << "] as " << rle
           << " top-left (" << x0 << "," << y0 << ")";
    }
    return ss.str();
  }

  // ---- escaped spaceships (spaceship-cutoff) ----
  // Index just after the longest circular run of unset bits (the "left edge" of a
  // compact set of rows/columns on the torus), or -1 if all/none are set.
  static int CircularLeftEdge(row_t occupied) {
    if (occupied == 0 || occupied == ROW_ALL)
      return -1;
    int start = RowCtz(occupied), best = -1, bestLen = -1, run = 0;
    for (int k = 1; k <= N; k++) {
      int i = (start + k) % N;
      if (!((occupied >> i) & 1)) {
        run++;
      } else {
        if (run > bestLen) { bestLen = run; best = i; }
        run = 0;
      }
    }
    return best;
  }

  // Moves a small object so its bounding box starts at (0, 0); same shape at a
  // different place gives the same result. False if it is not compact.
  static bool Normalise(LifeState &st) {
    row_t rows = 0;
    for (int i = 0; i < N; i++)
      rows |= st.state[i];
    int cx = CircularLeftEdge(NonEmptyColumns(st.state)), cy = CircularLeftEdge(rows);
    if (cx < 0 || cy < 0)
      return false;
    st.Move(-cx, -cy);
    st.gen = 0;
    return true;
  }

  // True if `obj`, evolved on its own, comes back as a moved copy of itself
  // within 32 generations. Results are cached per thread by shape.
  static bool IsSpaceship(const LifeState &obj) {
    LifeState key = obj;
    if (!Normalise(key))
      return false;
    static thread_local std::unordered_map<uint64_t, std::vector<std::pair<LifeState, bool>>> cache;
    uint64_t h = key.GetHash();
    auto &bucket = cache[h];
    for (const auto &e : bucket)
      if (e.first == key)
        return e.second;

    bool ship = false;
    const unsigned pop = obj.GetPop();
    LifeState cur = obj;
    for (int p = 1; p <= 32; p++) {
      cur.Step();
      unsigned cp = cur.GetPop();
      if (cp == 0 || cp > 4 * pop + 16)
        break;
      if (cp != pop)
        continue;
      LifeState n = cur;
      if (Normalise(n) && n == key) {
        ship = !(cur == obj); // same place: an oscillator, not a spaceship
        break;
      }
    }
    if (cache.size() > 100000)
      cache.clear();
    cache[h].emplace_back(key, ship);
    return ship;
  }

  // Cells within Chebyshev distance d of `st` (a (2d+1)x(2d+1) square around
  // each cell), in O(log d) shift rounds.
  static LifeState Dilate(const LifeState &st, int d) {
    LifeState a = st, b;
    for (int r = 0; r < d;) {
      int sh = std::min(2 * r + 1, d) - r; // radius r -> r + sh
      for (int i = 0; i < N; i++)
        b.state[i] = a.state[i] | RotateLeft(a.state[i], sh) | RotateRight(a.state[i], sh);
      for (int i = 0; i < N; i++)
        a.state[i] = b.state[i] | b.state[(i + sh) % N] | b.state[(i + N - sh) % N];
      r += sh;
    }
    a.RecalculateMinMax();
    return a;
  }

  // Cells of `active` (the non-catalyst cells) that belong to spaceships more
  // than `cutoff` cells away from all other active cells, provided something
  // other than those spaceships is left (the reaction or its debris).
  // `changing` = cells that differ in the next generation; only clusters that
  // contain some are looked at, since a still life is never a spaceship.
  static LifeState FarSpaceships(const LifeState &active, const LifeState &changing, int cutoff) {
    const unsigned maxShipPop = 150; // bigger clusters are taken to be the reaction
    LifeState far, seeds = changing & active;
    for (int tries = 0; tries < 8 && !seeds.IsEmpty(); tries++) {
      std::pair<int, int> c = seeds.FirstOn();
      LifeState comp;
      comp.Set(c.first, c.second);
      bool small = true;
      while (true) {
        LifeState grown = Dilate(comp, cutoff) & active;
        if (grown == comp)
          break;
        comp = grown;
        if (comp.GetPop() > maxShipPop) {
          small = false;
          break;
        }
      }
      seeds &= ~(small ? comp : Dilate(comp, cutoff));
      if (small && IsSpaceship(comp))
        far |= comp;
    }
    if (far.IsEmpty() || (active & ~far).IsEmpty())
      return LifeState();
    return far;
  }

  // Union of all search-area rectangles (default rectangle if none were
  // given) minus all exclude-area rectangles.
  LifeState BuildSearchBounds(bool print = true) const {
    LifeState bounds;
    if (params.searchAreas.empty() && !params.catalystZones.empty()) {
      bounds = ~LifeState(); // catalyst-zone limits the placements instead
    } else if (params.searchAreas.empty()) {
      bounds = LifeState::SolidRect(params.searchArea[0], params.searchArea[1],
                                    params.searchArea[2], params.searchArea[3]);
    } else {
      for (const auto &r : params.searchAreas)
        bounds |= LifeState::SolidRect(r[0], r[1], r[2], r[3]);
    }
    for (const auto &r : params.excludeAreas)
      bounds &= ~LifeState::SolidRect(r[0], r[1], r[2], r[3]);
    bounds &= FundamentalDomain(params.symmetry);

    if (print)
      std::cout << "Search area: " << std::max<size_t>(1, params.searchAreas.size())
              << " rectangle(s), " << params.excludeAreas.size() << " excluded, "
              << bounds.GetPop() << " positions" << std::endl;
    return bounds;
  }

  // --merge: reads the .cfdata files of earlier runs of this same input file,
  // combines them (duplicates removed as in a search) and writes the usual
  // output files with "_merged" added to their names.
  int Merge(const std::vector<std::string> &files) {
    for (std::string *f : {&params.outputFile, &params.fullReportFile, &params.outputList, &params.uniqueFile})
      if (!f->empty())
        *f = AddSuffix(*f, "_merged");
    params.fullReportFile.clear(); // the data files hold the filtered results only
    LifeState patSym;
    patSym.JoinWSymChain(pat, params.symmetryChain);
    unsigned read = 0, files_ok = 0;
    bool allFinished = true;
    for (const std::string &fn : files) {
      std::ifstream in(fn.c_str());
      if (!in.good())
        Die("--merge: cannot open " + fn);
      std::string line;
      bool finished = false, headerOk = false;
      while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r')
          line.pop_back();
        if (line.compare(0, 8, "# input ") == 0) {
          headerOk = true;
          finished = line.find("| finished yes") != std::string::npos;
          std::string want = "| rule " + RuleForHeader() + " | torus " + std::to_string(N) + " ";
          if (line.find(want) == std::string::npos)
            std::cout << "Warning: " << fn << " was made with a different rule or torus (" << line.substr(2) << ")" << std::endl;
          continue;
        }
        if (line.empty() || line[0] == '#')
          continue;
        std::vector<std::string> f;
        split(line, '\t', f);
        if (f.size() < 7)
          Die("--merge: " + fn + " has a damaged line: " + line.substr(0, 60));
        std::vector<std::string> c = splitwhitespace(f[1]), gl = splitwhitespace(f[2]);
        if (c.size() != 3 || gl.size() != 3)
          Die("--merge: " + fn + " has a damaged line: " + line.substr(0, 60));
        SolutionInfo info;
        info.cats = LifeState::Parse(c[0].c_str(), atoi(c[1].c_str()), atoi(c[2].c_str()));
        LifeState goal;
        if (gl[0] != "-")
          goal = LifeState::Parse(gl[0].c_str(), atoi(gl[1].c_str()), atoi(gl[2].c_str()));
        if (f[3] != "-") {
          std::vector<std::string> u;
          split(f[3], ',', u);
          for (const std::string &x : u)
            info.usage.push_back(atoi(x.c_str()));
        }
        info.catCount = atoi(f[4].c_str());
        info.edge = f[5] == "1";
        info.desc = f[6];
        for (size_t k = 7; k < f.size(); k++)
          info.desc += "\t" + f[k];
        const unsigned firstGen = atoi(f[0].c_str());
        LifeState init = info.cats | patSym;
        LifeState after = init, stepped = info.cats;
        after.Step(firstGen);
        stepped.Step(firstGen);
        Configuration conf;
        AddSolution(init, after, stepped, conf, firstGen, goal, nullptr, info);
        read++;
      }
      if (!headerOk)
        Die("--merge: " + fn + " is not a CatForce .cfdata file");
      allFinished = allFinished && finished;
      files_ok++;
    }
    searchFinished = true;
    FinalReport();
    std::cout << "Merged " << files_ok << " file(s): " << read << " solutions read, " << found << " kept"
              << (params.keepDuplicates ? "" : " (" + std::to_string(categoryContainer->duplicatesDropped) + " duplicates removed)")
              << (allFinished ? "" : "; note: some of these runs had not finished") << std::endl;
    return 0;
  }

  // countOnly: run only the first level and return how many first placements
  // this run (part / resume-from) will search.
  unsigned long long Search(bool countOnlyPass = false) {
    countOnly = countOnlyPass;
    countedRoots = 0;
    Configuration config;
    config.count = 0;
    config.transparentCount = 0;
    config.mustIncludeCount = 0;
    config.state.JoinWSymChain(pat, params.symmetryChain);

    LifeState bounds = BuildSearchBounds(!countOnlyPass); // already limited to the symmetry's fundamental domain

    if (bounds.IsEmpty()) {
      std::cout << "The search area is empty (check search-area / exclude-area / symmetry)" << std::endl;
      exit(1);
    }

    // catalyst-zone / no-catalyst-zone: forbid origins where any cell of the
    // catalyst (in any phase) would lie outside / inside them
    LifeState zone, noZone;
    zone.JoinWSymChain(catalystZone, params.symmetryChain);
    noZone.JoinWSymChain(noCatalystZone, params.symmetryChain);
    std::vector<LifeState> masks(catalysts.size());
    for (unsigned s = 0; s < catalysts.size(); s++) {
      masks[s] = config.state.Convolve(catalysts[s].reactionMask) | ~bounds;
      if (!catalystZone.IsEmpty() || !noCatalystZone.IsEmpty() || !catalysts[s].zone.IsEmpty()) {
        LifeState cells;
        for (const LifeState &ph : catalysts[s].phases)
          cells |= ph;
        cells.Transform(Rotate180OddBoth);
        cells.RecalculateMinMax();
        if (!catalystZone.IsEmpty())
          masks[s] |= (~zone).Convolve(cells);
        if (!noCatalystZone.IsEmpty())
          masks[s] |= noZone.Convolve(cells);
        if (!catalysts[s].zone.IsEmpty()) { // cat option zone
          LifeState own;
          own.JoinWSymChain(catalysts[s].zone, params.symmetryChain);
          masks[s] |= (~own).Convolve(cells);
        }
      }
    }

    // fixed-cat: in place from generation 0, like catalysts the search placed
    LifeState required = alsoRequired, antirequired;
    antirequired.JoinWSymChain(keepEmpty, params.symmetryChain);
    std::vector<TargetSet> baseTargets(totalCats);
    for (unsigned k = 0; k < fixedIdx.size(); k++) {
      const CatalystData &c = catalysts[fixedIdx[k]];
      const int x = fixedPos[k].first, y = fixedPos[k].second;
      LifeState placed = c.state;
      placed.Move(x, y);
      LifeState sym;
      sym.JoinWSymChain(placed, params.symmetryChain);
      if (!(config.state & sym).IsEmpty())
        Die("fixed-cat " + c.name + " at (" + std::to_string(x) + "," + std::to_string(y) + ") overlaps the pattern or another fixed-cat");
      config.curx[config.count] = x;
      config.cury[config.count] = y;
      config.curs[config.count] = fixedIdx[k];
      config.count++;
      config.startingCatalysts |= sym;
      config.catalystsState |= sym;
      config.state |= sym;
      required.Join(c.required, x, y);
      antirequired.Join(c.antirequired, x, y);
      if (!c.isBlinker) {
        LifeTarget t;
        t.wanted = placed;
        t.unwanted = c.target.unwanted;
        t.unwanted.Move(x, y);
        baseTargets[k].Add(t);
      } else {
        for (const LifeState &phase : c.phases) {
          LifeState ph = phase;
          LifeTarget t(ph);
          t.wanted.Move(x, y);
          t.unwanted.Move(x, y);
          baseTargets[k].Add(t);
        }
      }
      for (unsigned t = 0; t < catalysts.size(); t++)
        masks[t].Join(catalystCollisionMasks[fixedIdx[k] * catalysts.size() + t], x, y);
    }
    {
      LifeState next = config.catalystsState;
      next.Step();
      config.catalystsStill = next == config.catalystsState;
    }

    unsigned threads = params.numThreads;
    if (threads == 0)
      threads = std::max(1u, std::thread::hardware_concurrency());
    if (countOnlyPass)
      threads = 1;

    gNumWorkers = threads;
    gCurRoot.reset(new std::atomic<unsigned long>[threads]);
    for (unsigned t = 0; t < threads; t++) gCurRoot[t] = 0;
    auto worker = [&](unsigned part) {
      tlsPart = part;
      tlsParts = threads;
      tlsRootCounter = 0; tlsTicket = ~0UL; tlsQual = 0;
      std::vector<LifeState> myMasks = masks;
      std::vector<TargetSet> shiftedTargets = baseTargets;
      Configuration myConfig = config;
      RecursiveSearch(myConfig, myConfig.state, required, antirequired, std::move(myMasks), shiftedTargets,
                      std::array<unsigned, MAX_CATALYSTS>(), std::array<unsigned, MAX_CATALYSTS>(), nullptr, 0);
      if (!stopFlag) gCurRoot[part] = ~0UL; // this thread finished all its share
    };

    if (threads == 1) {
      worker(0);
    } else {
      std::cout << "Using " << threads << " threads" << std::endl;
      std::vector<std::thread> pool;
      for (unsigned t = 0; t < threads; t++)
        pool.emplace_back(worker, t);
      for (auto &th : pool)
        th.join();
    }
    countOnly = false;
    if (countOnlyPass) {
      interactingRejected = 0; // counted again by the real search
      progressGen = 0;
    }
    return countedRoots.load();
  }

  // Printed by --check and at the start of a search.
  void PrintSummary() const {
    int x0, y0, w, h;
    BoundingBoxRle(pat, x0, y0, w, h);
    std::cout << "Summary" << std::endl
              << "  rule:          " << RuleDescription() << std::endl
              << "  torus:         " << N << "x" << N << std::endl
              << "  spaceships:    " << (params.spaceshipCutoff > 0 ? "no catalysts for spaceships more than " + std::to_string(params.spaceshipCutoff) + " cells from the reaction" : "no cutoff") << std::endl
              << "  pattern:       " << pat.GetPop() << " cells, " << w << "x" << h << " box at (" << x0 << "," << y0 << ")" << std::endl
              << "  symmetry:      " << params.symmetryName << std::endl
              << "  baits:         " << params.baits.size() << std::endl
              << "  fixed-cat:     " << params.fixedCats.size() << std::endl
              << "  zones:         " << (catalystZone.IsEmpty() ? std::string("-") : std::to_string(catalystZone.GetPop()) + " catalyst-zone cells")
              << (noCatalystZone.IsEmpty() ? "" : ", " + std::to_string(noCatalystZone.GetPop()) + " no-catalyst-zone cells")
              << (keepEmpty.IsEmpty() ? "" : ", " + std::to_string(keepEmpty.GetPop()) + " keep-empty cells") << std::endl
              << "  catalysts:     " << inputNames.size() << " cat line(s), " << catalysts.size()
              << " placements to try (orientations/phases included)" << std::endl
              << "  to place:      " << (params.minCatalysts < params.numCatalysts ? std::to_string(params.minCatalysts) + " to " : std::string())
              << params.numCatalysts << " catalyst(s), stable for " << params.stableInterval << " gens" << std::endl
              << "  generations:   max " << params.maxGen << ", first interaction " << params.startGen << ".." << params.lastGen << std::endl
              << "  filters:       " << params.targetFilter.size() << std::endl
              << "  threads:       " << (params.numThreads == 0 ? std::string("auto (") + std::to_string(std::max(1u, std::thread::hardware_concurrency())) + ")" : std::to_string(params.numThreads)) << std::endl
              << "  output:        " << params.outputFile
              << (params.fullReportFile.empty() ? "" : ", " + params.fullReportFile + " (full report)")
              << (params.outputList.empty() ? "" : ", " + params.outputList + " (list)") << std::endl;
    BuildSearchBounds();
  }

  void PrintUsage() const {
    unsigned shown = 0;
    std::vector<std::pair<unsigned, size_t>> order;
    for (size_t i = 0; i < catalystUsage.size(); i++)
      if (catalystUsage[i] > 0)
        order.emplace_back(catalystUsage[i], i);
    if (order.empty())
      return;
    std::sort(order.begin(), order.end(), [](const std::pair<unsigned, size_t> &a, const std::pair<unsigned, size_t> &b) {
      return a.first != b.first ? a.first > b.first : a.second < b.second; });
    std::cout << "Catalyst usage (solutions using each cat line):" << std::endl;
    std::cout << std::setfill(' '); // PrintTime leaves '0' as the fill character
    for (auto &o : order) {
      if (shown++ >= 15) { std::cout << "  ..." << std::endl; break; }
      std::cout << "  " << std::setw(5) << o.first << "  line " << (o.second + 1) << "  " << inputNames[o.second] << std::endl;
    }
  }

  void
  RecursiveSearch(Configuration config, LifeState history, const LifeState &requiredMask, const LifeState &antirequired,
                  std::vector<LifeState> masks,
                  std::vector<TargetSet> &shiftedTargets, // This can be shared

                  std::array<unsigned, MAX_CATALYSTS> missingTime,
                  std::array<unsigned, MAX_CATALYSTS> recoveredTime,
                  Trajectory *parentTraj, row_t diffCols) {
    // This node's own evolution: stored if it can have children (they read it),
    // otherwise only stepped forward.
    const bool isLeaf = config.count == totalCats;
    Trajectory traj;
    traj.parent = parentTraj;
    traj.diff = diffCols;
    traj.start = config.state.gen;
    traj.limit = params.maxGen;
    if (!isLeaf) {
      traj.states.push_back(config.state);
      traj.changed.push_back(0);
    }
    auto advance = [&]() {
      if (isLeaf)
        traj.Advance(config.state);
      else
        config.state = traj.At(config.state.gen + 1);
    };
    bool success = false;
    bool failure = false;

    // Required/antirequired cells are rare; skip those tests when there are none.
    // (also-required counts: it is the initial requiredMask.)
    const bool checkRequired = !requiredMask.IsEmpty() || !antirequired.IsEmpty();
    unsigned successtime = 0;
    unsigned failuretime = 0;

    // Lightcone: simulate with no further placements to find the first failure
    // (time F, cells P). Any later catalyst must be placed at some gen g <= F
    // close enough to P to influence it by F (light speed: 1 cell/gen).
    int lcF = -1;
    LifeState lcP;
    if (params.lightcone && config.count > 0 && config.count < totalCats) {
      LifeState cs = config.catalystsState;
      auto miss = missingTime;
      for (unsigned g2 = config.state.gen; g2 < params.maxGen && lcF < 0; g2++, checkRequired ? cs.Step() : void()) {
        const LifeState &st = traj.At(g2);
        if (checkRequired) {
          if (!(((st ^ (cs | alsoRequired)) & requiredMask).IsEmpty()) || !st.AreDisjoint(antirequired)) {
            lcF = g2; lcP = requiredMask | antirequired; break;
          }
        }
        for (unsigned i = 0; i < config.count; i++) {
          if (st.gen % catalysts[config.curs[i]].period != 0) continue;
          if (shiftedTargets[i].MatchedBy(st)) miss[i] = 0;
          else if ((miss[i] += catalysts[config.curs[i]].period) > catalysts[config.curs[i]].maxDisappear) {
            lcF = g2;
            for (auto &t : shiftedTargets[i].alts)
              lcP |= params.lightconeTight ? ((t.wanted & ~st) | (t.unwanted & st)) : (t.wanted | t.unwanted);
          }
        }
      }
    }
    LifeState lcZone;
    int lcZoneGen = -1;

    for (unsigned g = config.state.gen; g < params.maxGen; g++) {
      if (ShouldStop())
        break;

      if (config.count == rootCount && g > params.lastGen)
        return;
      if (config.count == rootCount && tlsPart == 0)
        progressGen = g;

      if (config.count == rootCount && tlsPart == 0 && params.verbose) {
        std::cout << "Collision at gen " << g << std::endl;
      }

      if (checkRequired) {
        LifeState desiredState = config.catalystsState | alsoRequired;
        LifeState maskedDiff = (config.state ^ desiredState) & requiredMask;

        if (!maskedDiff.IsEmpty() || !config.state.AreDisjoint(antirequired)) {
          failure = true;
        }
      }

      for (unsigned i = 0; i < config.count; i++) {
        if (config.state.gen % catalysts[config.curs[i]].period != 0)
          continue;

        if (shiftedTargets[i].MatchedBy(config.state)) {
          missingTime[i] = 0;
          recoveredTime[i] += catalysts[config.curs[i]].period;
        } else {
          missingTime[i] += catalysts[config.curs[i]].period;
          recoveredTime[i] = 0;
        }

        if (missingTime[i] > catalysts[config.curs[i]].maxDisappear) {
          failure = true;
        }
      }

      // Child configurations re-process this same generation after placing a
      // catalyst, so they must inherit the bait counters from *before* it.
      const auto baitAbsentBefore = config.baitAbsent;
      const auto baitStreakBefore = config.baitStreak;
      for (unsigned b = 0; b < baitTargets.size(); b++) {
        if (baitTargets[b].MatchedBy(config.state)) {
          config.baitStreak[b]++;
        } else {
          config.baitStreak[b] = 0;
          config.baitAbsent[b]++;
        }
      }

      if (failure) {
        failuretime = g;
        break;
      }

      // Too late: a catalyst that is missing now needs ceil(stable-interval / period)
      // more checks in a row, and a bait that is missing needs `survive` more
      // generations, before generation max-gen - 1. If that no longer fits, no
      // result can come from here (or from catalysts placed later), so stop.
      if (!success) {
        bool tooLate = false;
        for (unsigned i = 0; i < config.count && !tooLate; i++) {
          const CatalystData &c = catalysts[config.curs[i]];
          if (c.sacrificial || missingTime[i] == 0 || g % c.period != 0)
            continue;
          const unsigned need = (params.stableInterval + c.period - 1) / c.period * c.period;
          tooLate = g + need > params.maxGen - 1;
        }
        for (unsigned b = 0; b < baitTargets.size() && !tooLate; b++)
          tooLate = config.baitStreak[b] == 0 && g + baitSurvive[b] > params.maxGen - 1;
        if (tooLate)
          break;
      }


      // Catalysts are still missing here but none can be placed in time any
      // more, so nothing below this point of the search can succeed.
      // (a configuration with enough catalysts for a result goes on, it just stops placing)
      if (config.count < minTotal && (int)g > latestPlacementGenAny)
        break;

      // Try adding a catalyst
      if (config.state.gen >= params.startGen && config.count != totalCats && (int)g <= latestPlacementGenAny) {
        const LifeState next = traj.At(config.state.gen + 1);

        LifeState activePart = (~history).ZOI() & config.state & ~config.catalystsState;

        // Escaped spaceships far from the reaction: no catalyst may be placed for
        // them, and the places they pass are blocked (a catalyst there would have
        // been hit by them).
        if (params.spaceshipCutoff > 0 && !activePart.IsEmpty()) {
          LifeState ships = FarSpaceships(config.state & ~config.catalystsState, config.state ^ next,
                                          params.spaceshipCutoff);
          if (!ships.IsEmpty()) {
            for (unsigned s = 0; s < catalysts.size(); s++)
              masks[s] |= ships.Convolve(catalysts[s].reactionMask | catalysts[s].phaseReactionMask[g % catalysts[s].period]);
            activePart &= ~ships;
          }
        }

          if(!activePart.IsEmpty()) {
          for (unsigned s = 0; s < catalysts.size(); s++) {
            if(catalysts[s].hasLocus) {
              LifeState hitLocations = activePart.Convolve(catalysts[s].locusAvoidMask);
              masks[s] |= hitLocations;
            }
          }

          if (lcF >= 0) {
            if (lcZoneGen != (int)g) {
              lcZone = lcP;
              for (int k = g; k < lcF; k++) lcZone = lcZone.ZOI();
              lcZoneGen = g;
            }
          }
          for (unsigned s = 0; s < catalysts.size(); s++) {
            if (catalysts[s].fixed)
              continue;
            if (catalysts[s].maxUses) {
              unsigned uses = 0;
              for (unsigned j = 0; j < config.count; j++)
                uses += catalysts[config.curs[j]].inputIndex == catalysts[s].inputIndex;
              if (uses >= catalysts[s].maxUses)
                continue;
            }
            if (config.transparentCount == params.numTransparent && catalysts[s].transparent)
              continue;
            if (config.count == totalCats - 1 && config.mustIncludeCount == 0 && !catalysts[s].mustInclude)
              continue;
            if ((int)g > latestPlacementGen[s])
              continue; // too late for this catalyst to recover before max-gen

            LifeState newPlacements = activePart.Convolve(catalysts[s].phaseReactionMask[g % catalysts[s].period]) & ~masks[s];

            while (!newPlacements.IsEmpty()) {
              // Do the placement
              auto newPlacement = newPlacements.FirstOn();

              LifeState shiftedCatalyst = catalysts[s].state;
              shiftedCatalyst.Move(newPlacement.first, newPlacement.second);

              LifeState symCatalyst;
              symCatalyst.JoinWSymChain(shiftedCatalyst, params.symmetryChain);
              const LifeState startingSym = symCatalyst; // phase 0, as the output shows it
              for (unsigned k = 0; k < g % catalysts[s].period; k++)
                symCatalyst.Step();
              if (lcF >= 0 && symCatalyst.AreDisjoint(lcZone)) {
                newPlacements.Erase(newPlacement.first, newPlacement.second);
                continue;
              }

              // Do a one-step lookahead to see if the catalyst interacts. Adding the
              // catalyst can only change the next generation within one column of
              // the catalyst, so only those columns are stepped (not the whole board).
              LifeState testState = config.state;
              testState |= symCatalyst;
              const row_t catCols = NonEmptyColumns(symCatalyst.state);
              const row_t nearCols = catCols | RotateLeft(catCols) | RotateRight(catCols);
              {
                row_t withCat[N], catAlone[N];
                testState.StepColumnsInto(nearCols, withCat);
                symCatalyst.StepColumnsInto(nearCols, catAlone);

                bool interacts = false;
                for (row_t bits = nearCols; bits; bits &= bits - 1) {
                  int i = RowCtz(bits);
                  if (((next.state[i] | catAlone[i]) ^ withCat[i]) != 0) {
                    interacts = true;
                    break;
                  }
                }

                if (!interacts) {
                  if (config.count == rootCount && tlsPart == 0 && params.verbose) {
                    std::cout << "Skipping catalyst " << s << " at "
                              << newPlacement.first << ", "
                              << newPlacement.second
                              << " (no interaction) "
                              << std::endl;
                  }

                  // Note: we deliberately don't set the mask,
                  // because it may turn out that a catalyst here
                  // interacts properly in a later generation.
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              Configuration newConfig = config;
              newConfig.baitAbsent = baitAbsentBefore;
              newConfig.baitStreak = baitStreakBefore;
              newConfig.count += 1;
              newConfig.curx[config.count] = newPlacement.first;
              newConfig.cury[config.count] = newPlacement.second;
              newConfig.curs[config.count] = s;
              if (catalysts[s].transparent)
                newConfig.transparentCount++;
              if (catalysts[s].mustInclude)
                newConfig.mustIncludeCount++;

              newConfig.startingCatalysts |= startingSym;
              newConfig.catalystsState |= symCatalyst;
              newConfig.state = testState;

              // Is the catalysts-only state still a fixed point? If it was before, only
              // the columns next to the new catalyst can have changed.
              // With a global symmetry, reject placements whose images interact with the
              // catalysts already there (or with each other). Position-only property, so
              // the placement is masked like the other permanent rejections.
              if (!params.symmetryChain.empty() && !params.keepInteracting) {
                bool coexist;
                {
                  unsigned P = CatalystPeriodLcm(catalysts, config);
                  P = P / gcd(P, catalysts[s].period) * catalysts[s].period;
                  coexist = CatalystsCoexist(newConfig.catalystsState, P);
                }
                if (!coexist) {
                  interactingRejected++;
                  masks[s].Set(newPlacement.first, newPlacement.second);
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              newConfig.catalystsStill = false;
              if (config.catalystsStill) {
                row_t out[N];
                newConfig.catalystsState.StepColumnsInto(nearCols, out);
                bool still = true;
                for (row_t bits = nearCols; bits; bits &= bits - 1) {
                  int i = RowCtz(bits);
                  if (out[i] != newConfig.catalystsState.state[i]) {
                    still = false;
                    break;
                  }
                }
                newConfig.catalystsStill = still;
              }

              LifeState newRequiredMask = requiredMask;
              newRequiredMask.Join(catalysts[s].required, newPlacement.first, newPlacement.second);

              LifeState newAntirequired = antirequired;
              newAntirequired.Join(catalysts[s].antirequired, newPlacement.first, newPlacement.second);

              // Skip the 4-step lookahead entirely when there is nothing it could reject
              // (no required/antirequired cells and no also-required pattern).
              if (!newAntirequired.IsEmpty() || !newRequiredMask.IsEmpty()) {
                LifeState lookahead = newConfig.state;
                lookahead.Step();
                lookahead.Step();
                lookahead.Step();
                lookahead.Step();

                bool requiredFailed = false;
                bool antirequiredFailed = false;
                if (!lookahead.AreDisjoint(newAntirequired)) {
                  antirequiredFailed = true;
                } else {
                  LifeState lookaheadcats = newConfig.catalystsState;
                  if (!newConfig.catalystsStill) {
                    lookaheadcats.Step();
                    lookaheadcats.Step();
                    lookaheadcats.Step();
                    lookaheadcats.Step();
                  }

                  LifeState desiredState = lookaheadcats | alsoRequired;
                  LifeState maskedDiff = (lookahead ^ desiredState) & newRequiredMask;

                  if (!maskedDiff.IsEmpty()) {
                    requiredFailed = true;
                  }
                }

                if (requiredFailed || antirequiredFailed) {
                  if (config.count == rootCount && tlsPart == 0 && params.verbose) {
                    std::cout << "Skipping catalyst " << s << " at "
                              << newPlacement.first << ", "
                              << newPlacement.second
                              << " (is destroyed) "
                              << std::endl;
                  }

                  masks[s].Set(newPlacement.first, newPlacement.second);
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              shiftedTargets[config.count].Clear();
              if (!catalysts[s].isBlinker) {
                LifeTarget t;
                t.wanted = shiftedCatalyst;
                t.unwanted = catalysts[s].target.unwanted;
                t.unwanted.Move(newPlacement.first, newPlacement.second);
                shiftedTargets[config.count].Add(t);
              } else {
                // Accept the catalyst in any of its phases
                for (const LifeState &phase : catalysts[s].phases) {
                  LifeState ph = phase;
                  LifeTarget t(ph);
                  t.wanted.Move(newPlacement.first, newPlacement.second);
                  t.unwanted.Move(newPlacement.first, newPlacement.second);
                  shiftedTargets[config.count].Add(t);
                }
              }

              bool soloFail = false;
              if (params.soloRecovery && !catalysts[s].checkRecovery && !catalysts[s].sacrificial) {
                // Contact filter: with no further help, the catalyst must come back at least once
                // within max-active gens. Bad first contacts (catalyst just explodes) fail this.
                LifeState la = newConfig.state;
                soloFail = true;
                for (unsigned k = 0; k <= catalysts[s].maxDisappear && soloFail; k++) {
                  la.Step();
                  if (shiftedTargets[config.count].MatchedBy(la)) soloFail = false;
                }
              }
              if (catalysts[s].checkRecovery || soloFail) {
                LifeState lookahead = newConfig.state;
                if (!soloFail) lookahead.Step(catalysts[s].maxDisappear);
                if (soloFail || !shiftedTargets[config.count].MatchedBy(lookahead)) {
                  if (config.count == rootCount && tlsPart == 0 && params.verbose) {
                    std::cout << "Skipping catalyst " << s << " at "
                              << newPlacement.first << ", "
                              << newPlacement.second
                              << " (failed to recover completely) "
                              << std::endl;
                  }

                  masks[s].Set(newPlacement.first, newPlacement.second);
                  newPlacements.Erase(newPlacement.first, newPlacement.second);
                  continue;
                }
              }

              if (config.count == rootCount && tlsPart == 0 && params.verbose) {
                std::cout << "Placing catalyst " << s << " at "
                          << newPlacement.first << ", " << newPlacement.second
                  //                                                              << ": " << newConfig.state.RLE()
                          << std::endl;
              }

              // If we just placed the last catalyst the child never places anything,
              // so it needs no masks at all (saves copying one LifeState per
              // catalyst type for most nodes of the search).
              std::vector<LifeState> newMasks;
              if (newConfig.count != totalCats) {
                newMasks = masks;
                LifeState bounds;
                if (params.maxW != -1) {
                  LifeState rect = LifeState::SolidRect(newPlacement.first - params.maxW,
                                                        newPlacement.second - params.maxH,
                                                        2 * params.maxW - 1,
                                                        2 * params.maxH - 1);
                  bounds.JoinWSymChain(rect, params.symmetryChain);

                  for (unsigned t = 0; t < catalysts.size(); t++) {
                    newMasks[t] |= ~bounds;
                  }
                }

                for (unsigned t = 0; t < catalysts.size(); t++) {
                  newMasks[t].Join(catalystCollisionMasks[s * catalysts.size() + t],
                                   newPlacement.first, newPlacement.second);
                }
              }

              // With several threads, each one explores only its share of the
              // top-level placements (every thread enumerates the same sequence).
              // With --part K/N, only every N-th of them belongs to this run at all.
              bool myShare = true;
              if (config.count == rootCount) {
                unsigned long r = tlsRootCounter++;
                myShare = r >= params.resumeFrom && (r % params.numParts) == params.partIndex;
                if (myShare && tlsParts > 1 && (params.numCatalysts == 1 || params.splitLevel < 2)) {
                  unsigned long idx = tlsQual++;   // index among the roots of this run (same in every thread)
                  if (tlsTicket == ~0UL) tlsTicket = gNextTicket.fetch_add(1);
                  myShare = (tlsTicket == idx);
                  if (myShare) tlsTicket = gNextTicket.fetch_add(1); // next ticket is always > idx
                }
                if (myShare && gCurRoot) gCurRoot[tlsPart] = r;
              } else if (config.count == rootCount + 1 && tlsParts > 1 && params.numCatalysts >= 2 && params.splitLevel >= 2) {
                // Several catalysts: every thread runs each first placement, and the threads
                // share out the second-level placements (much better balance with few roots).
                unsigned long idx = tlsQual++;
                if (tlsTicket == ~0UL) tlsTicket = gNextTicket.fetch_add(1);
                myShare = (tlsTicket == idx);
                if (myShare) tlsTicket = gNextTicket.fetch_add(1);
              }

              if (myShare && countOnly && config.count == rootCount) {
                countedRoots++; // counting pass: only count it
              } else if (myShare) {
                RecursiveSearch(newConfig, history, newRequiredMask, newAntirequired, std::move(newMasks),
                                shiftedTargets, missingTime, recoveredTime, &traj, catCols);
                if (config.count == rootCount) {
                  // (with split-level 2 every thread runs every first placement: count one thread's)
                  if (!(tlsParts > 1 && params.numCatalysts >= 2 && params.splitLevel >= 2) || tlsPart == 0)
                    rootPlacementsDone++;
                  MaybeReport();
                }
              }

              masks[s].Set(newPlacement.first, newPlacement.second);
              newPlacements.Erase(newPlacement.first, newPlacement.second);
            }
          }
          }
          history |= activePart;
          config.state = next;
      } else {
        // No need to update history, not used for anything
        // history |= config.state;
        advance();
      }
      // (a leaf only looks at the catalysts-only state for required/antirequired cells)
      if (config.catalystsStill || (isLeaf && !checkRequired))
        config.catalystsState.gen++; // fixed point, or not needed: stepping would change nothing
      else
        config.catalystsState.Step();

      // Still block the locations that are hit too early
      if (g < params.startGen && !masks.empty()) {
        for (unsigned s = 0; s < catalysts.size(); s++) {
          LifeState hitLocations = (config.state & ~config.catalystsState).Convolve(catalysts[s].phaseReactionMask[g % catalysts[s].period]);
          masks[s] |= hitLocations;
        }
      }

      if (config.count >= minTotal && !success && config.mustIncludeCount > 0) {
        bool allRecovered = true;
        for (unsigned i = 0; i < config.count; i++) {
          if (catalysts[config.curs[i]].sacrificial)
            continue;
          if (recoveredTime[i] < params.stableInterval || missingTime[i] > 0) {
            allRecovered = false;
          }
        }
        for (unsigned b = 0; b < baitTargets.size() && allRecovered; b++) {
          if (config.baitAbsent[b] < baitAbsentMin[b] || config.baitStreak[b] < baitSurvive[b])
            allRecovered = false;
        }
        if (allRecovered) {
          success = true;
          successtime = g;
        }
      }

      if (config.count == rootCount)
        MaybeReport();
    }

    if (!failure)
      failuretime = filterMaxGen;

    if (success)
      ReportSolution(config, successtime, failuretime);
  }
};


// ---- --compile-rule ---------------------------------------------------------
// Writes rule_kernels.h (next to the program) with a compiled kernel for the
// input file's rule plus up to 7 rules compiled before, rebuilds with make, and
// restarts the program with the same arguments. Returns only on failure.
static std::string DirOf(const std::string &path) {
  size_t slash = path.find_last_of('/');
  return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
}

static bool FileExists(const std::string &path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

static void CompileRuleAndRestart(int argc, char *argv[]) {
  std::string dir = DirOf(argv[0]);
  if (!FileExists(dir + "/Makefile") || !FileExists(dir + "/CatForce.cpp"))
    Die("--compile-rule: cannot find the CatForce sources (Makefile, CatForce.cpp) next to the program in \"" + dir +
        "\". Run the program from its source folder, e.g. ./CatForce your.in --compile-rule");

  std::vector<std::string> lines = ReadInputLines(argv[1]);
  ApplyRuleFromFile(lines);
  if (useFastLife) {
    std::cout << "B3/S23 already has its own fast stepper; nothing to compile." << std::endl;
    return;
  }
  const std::string newRule = currentRuleName;
  const std::string newKey = RuleTableKey();

  // rules compiled earlier (kept so switching between a few rules needs no rebuild)
  std::vector<std::string> rules{newRule};
  std::vector<std::string> keys{newKey};
  {
    std::ifstream old((dir + "/rule_kernels.h").c_str());
    std::string line;
    while (std::getline(old, line) && rules.size() < 8) {
      if (line.compare(0, 9, "// RULE: ") != 0)
        continue;
      std::string r = line.substr(9), err;
      if (!SetRule(r, err) || useFastLife)
        continue;
      std::string k = RuleTableKey();
      if (std::find(keys.begin(), keys.end(), k) != keys.end())
        continue;
      rules.push_back(r);
      keys.push_back(k);
    }
  }

  auto t0 = std::chrono::steady_clock::now();
  std::stringstream hdr, entries;
  hdr << "// Generated by CatForce --compile-rule: bit-parallel stepping kernels for the\n"
      << "// rules below. Delete this file (and run make) to use the table stepper only.\n";
  for (const std::string &r : rules)
    hdr << "// RULE: " << r << "\n";
  hdr << "#pragma once\n\n";
  for (size_t i = 0; i < rules.size(); i++) {
    std::string err;
    SetRule(rules[i], err);
    bool birth[256] = {false}, survive[256] = {false};
    for (int idx = 0; idx < 512; idx++) {
      int m = TableIndexToMask(idx);
      if ((idx >> 4) & 1)
        survive[m] = table[idx];
      else
        birth[m] = table[idx];
    }
    int ops = 0;
    std::string name = "RuleKernel" + std::to_string(i);
    hdr << rulecompiler::GenerateKernel(birth, survive, name, rules[i], ops) << "\n";
    entries << "  {\"" << keys[i] << "\", \"" << rules[i] << "\", &StepCompiledTmpl<" << name
            << "<RuleVec>>, &StepColsCompiledTmpl<" << name << "<RuleVec>>},\n";
    std::cout << "Kernel for " << rules[i].substr(0, 60) << ": " << ops << " operations" << std::endl;
  }
  hdr << "static const CompiledRuleEntry kCompiledRules[] = {\n" << entries.str() << "};\n";
  SaveTextAtomically(dir + "/rule_kernels.h", hdr.str());
  double genSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::cout << "Wrote " << dir << "/rule_kernels.h (" << rules.size() << " rule(s), " << std::fixed
            << std::setprecision(1) << genSec << " s). Rebuilding..." << std::endl;

  const char *makeCmd = getenv("CATFORCE_MAKE");
  std::string cmd = std::string(makeCmd ? makeCmd : "make") + " -C '" + dir + "'";
  auto t1 = std::chrono::steady_clock::now();
  int status = system(cmd.c_str());
  double buildSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
  if (status != 0)
    Die("--compile-rule: the rebuild failed (\"" + cmd + "\"). rule_kernels.h was written; fix the build "
        "and run make yourself, or delete rule_kernels.h to go back.");
  std::cout << "Rebuilt in " << std::setprecision(1) << buildSec << " s. Restarting with the compiled kernel..."
            << std::endl
            << std::endl;

  std::vector<char *> args;
  for (int i = 0; i < argc; i++)
    if (std::string(argv[i]) != "--compile-rule")
      args.push_back(argv[i]);
  args.push_back(nullptr);
  execv(argv[0], args.data());
  Die("--compile-rule: rebuilt, but could not restart the program; just run it again.");
}

// "torus <size>" in the input file: if this binary was built for another size,
// build (if needed) and run the one for that size (CatForce / CatForce128)
// with the same arguments. Returns if this binary already fits.
static void SwitchTorusIfNeeded(int argc, char *argv[]) {
  int want = N;
  for (const std::string &line : ReadInputLines(argv[1])) {
    std::vector<std::string> elems = splitwhitespace(line);
    if (elems.size() >= 2 && elems[0] == "torus")
      want = ParseIntOrDie(elems[1], "torus", 1);
  }
  if (want == N)
    return;
  if (want != 64 && want != 128)
    Die("torus: " + std::to_string(want) + " is not supported (use 64 or 128)");
  const std::string dir = DirOf(argv[0]);
  const std::string name = want == 64 ? "CatForce" : "CatForce" + std::to_string(want);
  const std::string exe = dir + "/" + name;
  if (FileExists(dir + "/Makefile") && FileExists(dir + "/CatForce.cpp")) {
    if (!FileExists(exe))
      std::cout << "torus " << want << ": building " << exe << " (first use only)..." << std::endl;
    const char *makeCmd = getenv("CATFORCE_MAKE");
    std::string cmd = std::string(makeCmd ? makeCmd : "make") + " -s -C '" + dir + "' " + name;
    if (system(cmd.c_str()) != 0 && !FileExists(exe))
      Die("torus " + std::to_string(want) + ": building " + name + " failed (\"" + cmd + "\")");
  } else if (!FileExists(exe)) {
    Die("torus " + std::to_string(want) + " needs " + name + " next to this program (build it with: make " + name + ")");
  }
  std::cout << "torus " << want << ": running " << exe << std::endl;
  std::vector<char *> args(argv, argv + argc);
  std::string exeArg = exe;
  args[0] = &exeArg[0];
  args.push_back(nullptr);
  execv(exe.c_str(), args.data());
  Die("torus " + std::to_string(want) + ": could not start " + exe);
}

int main(int argc, char *argv[]) {
  if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
    std::cout << "Usage: CatForce <in file> [--check] [--threads N] [--max-time N[s|m|h]] [--max-results N] [--part K/N] [--verbose]" << std::endl
              << "  --check        parse the file, print a summary and exit" << std::endl
              << "  --show-input   print the keyword lines the file amounts to (layouts, << blocks expanded)" << std::endl
              << "  --merge F...   combine the .cfdata files of earlier runs of this file (parts, resumed runs)" << std::endl
              << "  --threads N    worker threads (0 = all cores)" << std::endl
              << "  --max-time T   stop after T (e.g. 90, 15m, 2h) and save" << std::endl
              << "  --max-results N  stop after N solutions" << std::endl
              << "  --part K/N     run only the K-th of N equal shares (e.g. one per machine)" << std::endl
              << "  --verbose      print every first-catalyst placement" << std::endl
              << "  --compile-rule build a fast kernel for this file's rule into the program (runs make,"
              << std::endl
              << "                 then starts the search); needed once per rule, up to 8 rules are kept" << std::endl
              << "  (the same options can be given in the file as threads / max-time / max-results)" << std::endl;
    exit(0);
  }

  bool checkOnly = false;
  int cliThreads = -1;
  long cliMaxTime = -1, cliMaxResults = -1;
  std::string cliPart;
  bool cliVerbose = false;
  bool compileRule = false;
  std::vector<std::string> mergeFiles;
  for (int i = 2; i < argc; i++) {
    std::string a = argv[i];
    auto next = [&](const char *what) -> std::string {
      if (i + 1 >= argc)
        Die(std::string(what) + " needs a value");
      return argv[++i];
    };
    if (a == "--check" || a == "--dry-run")
      checkOnly = true;
    else if (a == "--compile-rule")
      compileRule = true;
    else if (a == "--show-input") {
      for (const std::string &l : ReadInputLines(argv[1]))
        std::cout << l << std::endl;
      return 0;
    }
    else if (a == "--threads")
      cliThreads = ParseIntOrDie(next("--threads"), "--threads", 0);
    else if (a == "--max-time")
      cliMaxTime = ParseDuration(next("--max-time"), "--max-time");
    else if (a == "--max-results")
      cliMaxResults = ParseIntOrDie(next("--max-results"), "--max-results", 1);
    else if (a == "--part")
      cliPart = next("--part");
    else if (a == "--verbose")
      cliVerbose = true;
    else if (a == "--merge") {
      while (i + 1 < argc)
        mergeFiles.push_back(argv[++i]);
      if (mergeFiles.empty())
        Die("--merge needs one or more .cfdata files");
    }
    else
      Die("Unknown option: " + a);
  }

  if (compileRule) {
    CompileRuleAndRestart(argc, argv); // only returns when there is nothing to compile
  }
  SwitchTorusIfNeeded(argc, argv);

  std::cout << "Input: " << argv[1] << std::endl
            << "Initializing please wait..." << std::endl;

  CatalystSearcher searcher;
  searcher.skipMasks = !mergeFiles.empty();
  searcher.Init(argv[1]);
  if (!mergeFiles.empty())
    return searcher.Merge(mergeFiles);
  if (cliThreads >= 0) searcher.params.numThreads = cliThreads;
  if (cliMaxTime >= 0) searcher.params.maxTimeSec = (unsigned)cliMaxTime;
  if (cliMaxResults >= 0) searcher.params.maxResults = (unsigned)cliMaxResults;
  if (cliVerbose) searcher.params.verbose = true;
  if (!cliPart.empty())
    ParsePart(cliPart, searcher.params.partIndex, searcher.params.numParts);
  if (searcher.params.numParts > 1) {
    // separate files per part, so parts can run side by side in one directory
    std::string suffix = "_part" + std::to_string(searcher.params.partIndex + 1) + "of" +
                         std::to_string(searcher.params.numParts);
    searcher.params.outputFile = AddSuffix(searcher.params.outputFile, suffix);
    if (!searcher.params.fullReportFile.empty())
      searcher.params.fullReportFile = AddSuffix(searcher.params.fullReportFile, suffix);
    if (!searcher.params.outputList.empty())
      searcher.params.outputList = AddSuffix(searcher.params.outputList, suffix);
    if (!searcher.params.uniqueFile.empty())
      searcher.params.uniqueFile = AddSuffix(searcher.params.uniqueFile, suffix);
    std::cout << "Running part " << (searcher.params.partIndex + 1) << " of " << searcher.params.numParts
              << "; results go to " << searcher.params.outputFile << std::endl;
  }
  if (searcher.params.resumeFrom > 0) {
    std::string suffix = "_from" + std::to_string(searcher.params.resumeFrom);
    searcher.params.outputFile = AddSuffix(searcher.params.outputFile, suffix);
    if (!searcher.params.fullReportFile.empty()) searcher.params.fullReportFile = AddSuffix(searcher.params.fullReportFile, suffix);
    if (!searcher.params.outputList.empty()) searcher.params.outputList = AddSuffix(searcher.params.outputList, suffix);
    if (!searcher.params.uniqueFile.empty()) searcher.params.uniqueFile = AddSuffix(searcher.params.uniqueFile, suffix);
    std::cout << "Resuming after " << searcher.params.resumeFrom << " first placements; results go to "
              << searcher.params.outputFile << std::endl;
  }

  printf("Total elapsed time: %f seconds\n", searcher.ElapsedSeconds());

  if (checkOnly) {
    searcher.PrintSummary();
    std::cout << "  first placements: " << searcher.Search(true) << " (each is searched separately; the status line counts them)" << std::endl;
    std::cout << std::endl << "Check finished (no search was run)." << std::endl;
    return 0;
  }

  std::signal(SIGINT, OnInterrupt);
  std::cout << std::endl
            << "Initialization finished, searching..." << std::endl
            << std::endl;

  {
    double t0 = searcher.ElapsedSeconds();
    searcher.totalRoots = searcher.Search(true);
    searcher.countSeconds = searcher.ElapsedSeconds();
    std::cout << "First placements to search: " << searcher.totalRoots << " (counted in "
              << std::fixed << std::setprecision(1) << searcher.countSeconds - t0 << " s)" << std::endl;
  }
  searcher.Search();
  searcher.searchFinished = true;
  searcher.FinalReport();

  printf("\n\nFINISH\n");
  if (!searcher.stopReason.empty())
  {
    std::cout << "Stopped early: " << searcher.stopReason << " (results so far were saved)" << std::endl;
    unsigned long k = ~0UL;
    for (unsigned t = 0; t < gNumWorkers; t++) k = std::min(k, gCurRoot[t].load());
    if (k != ~0UL && k > searcher.params.resumeFrom) {
      std::cout << "To resume, add this line to the input file (new results go to files ending _from" << k
                << "):\n  resume-from " << k << std::endl;
      SaveTextAtomically(searcher.params.inputFile + ".resume", "resume-from " + std::to_string(k) + "\n");
    }
  }
  if (!searcher.params.keepDuplicates)
    printf("Duplicates removed: %u (filtered output), %u (full report)\n",
           searcher.categoryContainer->duplicatesDropped,
           searcher.fullCategoryContainer->duplicatesDropped);
  if (searcher.interactingRejected > 0)
    printf("Placements/solutions where catalysts (or their symmetric images) touched each other (%s): %u\n",
           searcher.params.keepInteracting ? "kept, keep-interacting is on" : "rejected",
           searcher.interactingRejected.load());
  if (searcher.edgeFlagged > 0)
    printf("Solutions reaching the torus edge (marked in output-list; \"torus-edge drop\" removes them, \"torus 128\" gives more room): %u\n",
           searcher.edgeFlagged.load());
  if (searcher.edgeDropped > 0)
    printf("Solutions dropped because they reach the torus edge: %u\n", searcher.edgeDropped.load());
  searcher.PrintUsage();
  printf("Total elapsed time: %f seconds\n", searcher.ElapsedSeconds());
}
