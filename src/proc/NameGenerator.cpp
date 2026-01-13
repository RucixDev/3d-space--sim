#include "stellar/proc/NameGenerator.h"

#include <array>
#include <cctype>
#include <string>

namespace stellar::proc {

/*
simplified control flow
no temp std:string objects 
a bit optimized
*/

static const std::array<const char*, 32> kSyllA = {
    "al","an","ar","as","be","ca","ce","da","de","el","en","er","es","fa","fi","ga",
    "ha","he","ia","in","is","ka","ke","la","le","li","ma","me","na","ne","or","ra"
};

static const std::array<const char*, 32> kSyllB = {
    "bar","bel","car","cer","dan","del","dor","far","fen","gar","gen","hal","hel","ian","iel","jor",
    "kal","kel","lar","len","mir","nal","ner","nor","par","pel","ran","rel","sar","sel","tor","ven"
};

static const std::array<const char*, 24> kSyllC = {
    "a","e","i","o","u","ae","ia","io","oa","ul","ur","on","en","is","as","os","ar","or","ir","ium","ara","ora","eus","is"
};

static std::string capitalize(std::string s) {
    if (!s.empty()) {
        s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    }
    return s;
}

std::string NameGenerator::systemName() {
    const int parts = rng_.range(2, 3);
    std::string s;

    if (parts == 3) {
        s += kSyllA[rng_.range<int>(0, static_cast<int>(kSyllA.size() - 1))];
    }
    s += kSyllB[rng_.range<int>(0, static_cast<int>(kSyllB.size() - 1))];
    s += kSyllC[rng_.range<int>(0, static_cast<int>(kSyllC.size() - 1))];

    return capitalize(s);
}

std::string NameGenerator::planetName(const std::string& sys, int index) {
    static const char* romans[] = {
        "I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII"
    };
    const int r = (index >= 0 && index < 12) ? index : (index % 12);
    return sys + " " + romans[r];
}

std::string NameGenerator::stationName(const std::string& sys, int index) {
    static const char* greek[] = {
        "Alpha","Beta","Gamma","Delta","Epsilon","Zeta","Eta","Theta","Iota","Kappa",
        "Lambda","Mu","Nu","Xi","Omicron","Pi","Rho","Sigma","Tau","Upsilon",
        "Phi","Chi","Psi","Omega"
    };
    const int g = (index >= 0) ? (index % static_cast<int>(sizeof(greek) / sizeof(greek[0]))) : 0;
    return sys + " Station " + greek[g];
}

} // namespace stellar::proc
