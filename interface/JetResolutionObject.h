// JetResolutionObject.h
#ifndef JETRESOLUTIONOBJECT_H
#define JETRESOLUTIONOBJECT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <initializer_list>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>   // ← Add this line
#include <stdexcept>
#include <TFormula.h>

template <typename T>
T clip(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
}

namespace JME {

enum class Variation {
    NOMINAL = 0,
    DOWN = 1,
    UP = 2
};

template <typename T, typename U>
struct bimap {
    std::unordered_map<T, U> left;
    std::unordered_map<U, T> right;

    bimap(std::initializer_list<typename std::unordered_map<T, U>::value_type> l) {
        for (auto& v : l) {
            left.insert(v);
            right.insert({v.second, v.first});
        }
    }

    bimap() = default;
    bimap(bimap&& rhs) noexcept {
        left = std::move(rhs.left);
        right = std::move(rhs.right);
    }
};

enum class Binning {
    JetPt = 0,
    JetEta,
    JetAbsEta,
    JetE,
    JetArea,
    Mu,
    Rho,
    NPV,
};

} // namespace JME

namespace std {
template<>
struct hash<JME::Binning> {
    size_t operator()(const JME::Binning& b) const noexcept {
        return hash<int>()(static_cast<int>(b));
    }
};
}

namespace JME {

class JetParameters {
public:
    using value_type = std::unordered_map<Binning, float>;

    JetParameters() = default;
    JetParameters(JetParameters&& rhs) noexcept;
    JetParameters(std::initializer_list<value_type::value_type> init);

    JetParameters& set(const Binning& bin, float value);

    static const bimap<Binning, std::string> binning_to_string;

    std::vector<float> createVector(const std::vector<Binning>& binning) const;

    JetParameters& setJetPt(float pt) {
        m_values[Binning::JetPt] = pt;
        return *this;
    }

    JetParameters& setJetEta(float eta) {
        m_values[Binning::JetEta] = eta;
        return *this;
    }

    JetParameters& setRho(float rho) {
        m_values[Binning::Rho] = rho;
        return *this;
    }

private:
    value_type m_values;
};

class JetResolutionObject {
public:
    struct Range {
        float min;
        float max;
        Range() = default;
        Range(float min_, float max_) : min(min_), max(max_) {}
        bool is_inside(float value) const { return value >= min && value <= max; }
    };

    class Definition {
    public:
        Definition() = default;
        explicit Definition(const std::string& definition);

        const std::vector<std::string>& getBinsName() const { return m_bins_name; }
        const std::vector<std::string>& getVariablesName() const { return m_variables_name; }
        const std::vector<Binning>& getBins() const { return m_bins; }
        const std::vector<Binning>& getVariables() const { return m_variables; }
        size_t nBins() const { return m_bins_name.size(); }
        size_t nVariables() const { return m_variables.size(); }

        void init();

        TFormula* getFormula() const { return m_formula.get(); }

    private:
        std::vector<std::string> m_bins_name;
        std::vector<std::string> m_variables_name;
        std::shared_ptr<TFormula> m_formula;
        std::vector<Binning> m_bins;
        std::vector<Binning> m_variables;
        std::string m_formula_str;
    };

    class Record {
    public:
        Record() = default;
        explicit Record(const std::string& record, const Definition& def);

        const std::vector<Range>& getBinsRange() const { return m_bins_range; }
        const std::vector<Range>& getVariablesRange() const { return m_variables_range; }
        const std::vector<float>& getParametersValues() const { return m_parameters_values; }

    private:
        std::vector<Range> m_bins_range;
        std::vector<Range> m_variables_range;
        std::vector<float> m_parameters_values;
    };

    JetResolutionObject();
    explicit JetResolutionObject(const std::string& filename);
    JetResolutionObject(const JetResolutionObject& other);
    JetResolutionObject& operator=(const JetResolutionObject& other);
    void dump() const;

    const Record* getRecord(const JetParameters& bins) const;
    float evaluateFormula(const Record& record, const JetParameters& variables) const;

private:
    Definition m_definition;
    std::vector<Record> m_records;
    bool m_valid = false;
};

// === Inline Implementations ===

inline JetParameters::JetParameters(JetParameters&& rhs) noexcept {
    m_values = std::move(rhs.m_values);
}

inline JetParameters::JetParameters(std::initializer_list<value_type::value_type> init) {
    for (const auto& i : init) {
        m_values.emplace(i);
    }
}

inline JetParameters& JetParameters::set(const Binning& bin, float value) {
    m_values[bin] = value;
    return *this;
}

inline std::vector<float> JetParameters::createVector(const std::vector<Binning>& binning) const {
    std::vector<float> values;
    for (const auto& bin : binning) {
        auto it = m_values.find(bin);
        if (it == m_values.end()) {
            throw std::runtime_error("Missing parameter for binning variable.");
        }
        values.push_back(it->second);
    }
    return values;
}


inline JetResolutionObject::JetResolutionObject() = default;

inline JetResolutionObject::JetResolutionObject(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Cannot open JER file: " + filename);
    }

    for (std::string line; std::getline(infile, line); ) {
        if (line.empty() || line[0] == '#') continue;

        size_t first = line.find('{');
        size_t last = line.find('}');

        if (first != std::string::npos && last != std::string::npos && first < last) {
            std::string def = line.substr(first + 1, last - first - 1);
            m_definition = Definition(def);
        } else {
            m_records.emplace_back(line, m_definition);
        }
    }

    m_valid = true;
}

inline JetResolutionObject::JetResolutionObject(const JetResolutionObject& other)
    : m_definition(other.m_definition),
      m_records(other.m_records),
      m_valid(other.m_valid) {}

inline JetResolutionObject& JetResolutionObject::operator=(const JetResolutionObject& other) {
    if (this != &other) {
        m_definition = other.m_definition;
        m_records = other.m_records;
        m_valid = other.m_valid;
    }
    return *this;
}

inline void JetResolutionObject::dump() const {
    std::cout << "JetResolutionObject with " << m_records.size() << " records.\n";
}

inline const JetResolutionObject::Record* JetResolutionObject::getRecord(const JetParameters& bins) const {
    if (!m_valid) return nullptr;
    auto bin_values = bins.createVector(m_definition.getBins());
    for (const auto& record : m_records) {
        bool match = true;
        const auto& ranges = record.getBinsRange();
        for (size_t i = 0; i < ranges.size(); ++i) {
            if (!ranges[i].is_inside(bin_values[i])) {
                match = false;
                break;
            }
        }
        if (match) return &record;
    }
    return nullptr;
}

inline float JetResolutionObject::evaluateFormula(const Record& record, const JetParameters& variables) const {
    if (!m_valid) return 1.f;
    auto formula = m_definition.getFormula();
    if (!formula) return 1.f;
    auto variable_values = variables.createVector(m_definition.getVariables());
    for (size_t i = 0; i < record.getParametersValues().size(); ++i) {
        formula->SetParameter(i, record.getParametersValues()[i]);
    }
    double vars[4] = {0};
    for (size_t i = 0; i < variable_values.size(); ++i) {
        vars[i] = clip(variable_values[i], record.getVariablesRange()[i].min, record.getVariablesRange()[i].max);
    }
    return formula->EvalPar(vars);
}

inline JetResolutionObject::Definition::Definition(const std::string& definition) {
    std::istringstream iss(definition);
    size_t nBins, nVars;
    iss >> nBins;
    m_bins_name.resize(nBins);
    for (size_t i = 0; i < nBins; ++i)
        iss >> m_bins_name[i];

    iss >> nVars;
    m_variables_name.resize(nVars);
    for (size_t i = 0; i < nVars; ++i)
        iss >> m_variables_name[i];

    iss >> m_formula_str;
    if (m_formula_str == "None" || m_formula_str == "none")
        m_formula_str.clear();

    init();
}

inline void JetResolutionObject::Definition::init() {
    if (!m_formula_str.empty())
        m_formula = std::make_shared<TFormula>("formula", m_formula_str.c_str());

    for (const auto& name : m_bins_name) {
        if (name == "JetPt") m_bins.push_back(Binning::JetPt);
        else if (name == "JetEta") m_bins.push_back(Binning::JetEta);
        else if (name == "JetAbsEta") m_bins.push_back(Binning::JetAbsEta);
        else if (name == "JetE") m_bins.push_back(Binning::JetE);
        else if (name == "JetArea") m_bins.push_back(Binning::JetArea);
        else if (name == "Mu") m_bins.push_back(Binning::Mu);
        else if (name == "Rho") m_bins.push_back(Binning::Rho);
        else if (name == "NPV") m_bins.push_back(Binning::NPV);
        else throw std::runtime_error("Unknown binning variable: " + name);
    }

    for (const auto& name : m_variables_name) {
        if (name == "JetPt") m_variables.push_back(Binning::JetPt);
        else if (name == "JetEta") m_variables.push_back(Binning::JetEta);
        else if (name == "JetAbsEta") m_variables.push_back(Binning::JetAbsEta);
        else if (name == "JetE") m_variables.push_back(Binning::JetE);
        else if (name == "JetArea") m_variables.push_back(Binning::JetArea);
        else if (name == "Mu") m_variables.push_back(Binning::Mu);
        else if (name == "Rho") m_variables.push_back(Binning::Rho);
        else if (name == "NPV") m_variables.push_back(Binning::NPV);
        else throw std::runtime_error("Unknown variable: " + name);
    }
}

inline JetResolutionObject::Record::Record(const std::string& line, const Definition& def) {
    std::istringstream iss(line);
    for (size_t i = 0; i < def.nBins(); ++i) {
        float min, max;
        iss >> min >> max;
        m_bins_range.emplace_back(min, max);
    }
    size_t nParameters;
    iss >> nParameters;
    for (size_t i = 0; i < def.nVariables(); ++i) {
        float min, max;
        iss >> min >> max;
        m_variables_range.emplace_back(min, max);
        nParameters -= 2;
    }
    for (size_t i = 0; i < nParameters; ++i) {
        float param;
        iss >> param;
        m_parameters_values.push_back(param);
    }
}

} // namespace JME

#endif // JETRESOLUTIONOBJECT_H