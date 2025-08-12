#ifndef JETRESOLUTION_H
#define JETRESOLUTION_H

#include "JetResolutionObject.h"
#include <memory>
#include <string>
#include <vector>

namespace JME {

class JetResolution {
public:
    JetResolution() = default;

    explicit JetResolution(const std::string& filename)
        : m_object(std::make_shared<JetResolutionObject>(filename)) {}

    float getResolution(const JetParameters& parameters) const {
        const JetResolutionObject::Record* record = m_object->getRecord(parameters);
        if (!record) return 1.0f;
        return m_object->evaluateFormula(*record, parameters);
    }

private:
    std::shared_ptr<JetResolutionObject> m_object;
};

class JetResolutionScaleFactor {
public:
    JetResolutionScaleFactor() = default;

    explicit JetResolutionScaleFactor(const std::string& filename)
        : m_object(std::make_shared<JetResolutionObject>(filename)) {}

    float getScaleFactor(const JetParameters& parameters, Variation variation = Variation::NOMINAL) const {
        const JetResolutionObject::Record* record = m_object->getRecord(parameters);
        if (!record) return 1.0f;
        return record->getParametersValues().at(static_cast<size_t>(variation));
    }

private:
    std::shared_ptr<JetResolutionObject> m_object;
};

} // namespace JME

#endif // JETRESOLUTION_H