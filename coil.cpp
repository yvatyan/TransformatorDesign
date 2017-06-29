#include <cmath>
#include <stddef.h>

#include <iostream> //TODO: Remove after testing

#define SQR(x) ((x)*(x))

namespace constants {
    const double pi     = M_PI;
    const double u0     = 4*pi*1e-7;
    const double sqrt3  = sqrt(3.);
}
namespace aux {
    double KelvinToCelsius(double K) {
        return K + 273.15;
    }
    double CelsiusToKelvin(double C) {
        return C - 273.15;
    }
    template <typename ReturnType, typename ArgType>
    class Function {
        public:
            virtual ReturnType operator()(ArgType arg) = 0;
    };
    class FunctionElipticIntegralOfFirstKind : public Function<double, double> {
        private:
            double k;
        public:
            FunctionElipticIntegralOfFirstKind(double k) : k(k) {};
            double operator()(double arg) {
               return 1./sqrt(1-k*SQR(sin(arg)));
            }
    };   
    class FunctionElipticIntegralOfSecondKind : public Function<double, double> {
        private:
            double k;
        public:
            FunctionElipticIntegralOfSecondKind(double k) : k(k) {};
            double operator()(double arg) {
               return sqrt(1-k*SQR(sin(arg)));
            }
    };  
    template<typename ReturnType, typename ArgType>
    ReturnType CalculateDefineIntegralWithLeftRectangleMethod(
        Function<ReturnType, ArgType>* function,
        double lowerBorder,
        double upperBorder,
        size_t nodsNumber
    ) {
        double step = (upperBorder - lowerBorder)/nodsNumber;
        double Sum = 0.;
        for(size_t i = 0; i < nodsNumber-1; ++i) {
            Sum += (*function)(lowerBorder + i*step);
        }
        return static_cast<ReturnType>(Sum*step);
    }
    template<typename ReturnType, typename ArgType>
    ReturnType CalculateDefineIntegralWithSimpsonMethod(
        Function<ReturnType, ArgType>* function,
        double lowerBorder,
        double upperBorder,
        size_t nodsNumber
    ) {
        double step = (upperBorder - lowerBorder)/nodsNumber;
        double Sum = 0.;
        for(size_t i = 1; i < nodsNumber-1; ++i) {
            Sum += (*function)(lowerBorder + i*step)*(i&1 ? 4. : 2.);
        }
        return static_cast<ReturnType>((Sum + (*function)(lowerBorder) + (*function)(upperBorder))*step/3.);
    }
/****************************** Functions' Short names ******************************/
    double KtC(double K) {
    	return KelvinToCelsius(K);
    }
    double CtK(double C) {
    	return CelsiusToKelvin(C);
    }
/**********************************************************************************/
}
namespace CoilCalculation {
    class ElectricalEnvironment {
        public:
            double Current()       const;   // A
            double Voltage()       const;   // V
            double Frequancy()     const;   // Hz  
            double Temperature()   const;   // K
            void   SetCurrent(
	    		double i);
            void   SetVoltage(
	    		double v);
            void   SetFrequancy(
	    		double f);
            void   SetTemperature(
	    		double t);
        private:
            double current;
            double voltage;
            double frequancy;
            double temperature;
    };
    class Material {
        public:
            enum MaterialType {
	    	None = 0,
                Cu,
                Fe
            };
            Material();
            Material(MaterialType t);
            
            double       Density()         const;   // kg/m^3
            double       Resistivity()     const;   // Ohm*m [at 20C]
            double       HeatingCapacity() const;   // J/(kg*K)
            double       Permeability()    const;   // H/m
            double       MeltingPoint()    const;   // K
            bool         IsMetal()         const;
            MaterialType Type()            const;
        private:
            MaterialType type;
    };
    class Body {
        public:
            enum BodyType {
                Box,
                Cylinder
            };
            virtual double   BasePerimeter()           const = 0;    // m
            virtual double   BaseArea()                const = 0;    // m^2
            virtual double   SideArea()                const = 0;    // m^2
            virtual double   Height()                  const = 0;    // m
            virtual double   BaseIncircleRadius()      const = 0;    // m
            virtual double   BaseExcircleRadius()      const = 0;    // m
            virtual double   Volume()                  const = 0;    // m^3
            virtual double   BasePerimeterStretchedBy(
                                double amount)         const = 0;    // m
                    BodyType Type()                    const;
                    Material MaterialType()            const;
        protected:
            Body(BodyType t, Material m);
        private:
            BodyType type;
            Material material;
    };
    class BoxBody : public Body {
        public:
            BoxBody(double h, double l, double w, Material m);
            
            double   BasePerimeter()           const;
            double   BaseArea()                const;
            double   SideArea()                const;
            double   Height()                  const;
            double   Length()                  const; // m
            double   Width()                   const; // m
            double   BaseIncircleRadius()      const;
            double   BaseExcircleRadius()      const;
            double   Volume()                  const;
            double   BasePerimeterStretchedBy(
                        double eachSideAmount) const;
        private:
            double height;
            double length;
            double width;
    };
    class CylindricBody : public Body {
        public:
            CylindricBody(double h, double r, Material m);
            
            double   BasePerimeter()           const;
            double   BaseArea()                const;
            double   SideArea()                const;
            double   Height()                  const;
            double   Radius()                  const; // m
            double   BaseIncircleRadius()      const;
            double   BaseExcircleRadius()      const;
            double   Volume()                  const;
            double   BasePerimeterStretchedBy(
                        double radiusAmount)   const;
        private:
            double height;
            double radius;
    };
    class Wire {
       public:
            Wire(ElectricalEnvironment* e, double len, double metalD, double isolatedD, Material m);

            double   Length()                 		    const;  // m
            double   MetalDiameter()          		    const;  // m
            double   IsolatedDiameter()	      		    const;  // m
            double   CutArea(
	    		bool isolated = true) 		    const;  // m^2
            double   Mass(
	    		bool treatIsolationAsMetal = false) const;  // kg
            double   Resistance()             		    const;  // Ohm
            double   Heating()                		    const;  // J/s
            double   TimeToMelt()             		    const;  // s
            Material Metal()                  		    const;
        private:
            ElectricalEnvironment* env;
            double            length;
            double            metalDiameter;
            double            isolatedDiameter;
            Material          metal;
    };
    class Coil {
        public:
            enum PackingType {
                Quadratic,
                Compact
            };
            Coil(ElectricalEnvironment* e, Wire w, PackingType pack, CylindricBody core);
            Coil(ElectricalEnvironment* e, Wire w, PackingType pack, BoxBody core);
            ~Coil();

            size_t      Turns()            const;
            size_t      Layers()           const;
            double      WindingThickness() const;
            PackingType Packing()          const;
            Material    CoreMaterial()     const;
            double      Inductivity()      const;
            double      Resistance()       const;
        private:
            struct WindingParams {
                size_t turns;
                size_t layers;
                double thickness;
            }                params;
            ElectricalEnvironment* env;
            Wire              wire;
            Body*             coreBody; 
            PackingType       packing;
	    double	      inductivity;

            void winding();
	    double stretchAmount(
	    		double layer_iter) 	        const;
	    double turnInductivity(
	    		double distanceBetweenTurns,
			double firstTurnRadius,
			double secondTurnRadius = -1.);
    };
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace CoilCalculation {
/************************ class ElectricalEnvironment ************************/
double ElectricalEnvironment::Current() const {
    return current;    
}
double ElectricalEnvironment::Voltage() const {
    return voltage;
}
double ElectricalEnvironment::Frequancy() const {
    return frequancy;
}
double ElectricalEnvironment::Temperature() const {
    return temperature;
}
void ElectricalEnvironment::SetCurrent(double i) {
    current = i;
}
void ElectricalEnvironment::SetVoltage(double v) {
    voltage = v;
}
void ElectricalEnvironment::SetFrequancy(double f) {
    frequancy = f;
}
void ElectricalEnvironment::SetTemperature(double t) {
    temperature = t;
}
/**************************** class Material ****************************/
Material::Material() {
    type = None;
}
Material::Material(Material::MaterialType t) {
    type = t;
}
double Material::Density() const {
    switch(type) {
        case Cu :   return 8930.;
        case Fe :   return 7850.;
        default :   return -1;
    };
}
double Material::Resistivity() const {
    switch(type) {
        case Cu :   return 1.68e-8;
        case Fe :   return 9.71e-8;
        default :   return -1;
    };
}
double Material::HeatingCapacity() const {
    switch(type) {
        case Cu :   return 385;
        case Fe :   return 444;
        default :   return -1;
    };
}
double Material::Permeability() const {
    switch(type) {
        case Cu :   return 1.26e-6;
        case Fe :   return 6.30e-3;
        default :   return -1;
    };
}
double Material::MeltingPoint() const {
     switch(type) {
        case Cu :   return 1356.55;
        case Fe :   return 1812.15;
        default :   return false;
    };
}
bool Material::IsMetal() const {
     switch(type) {
        case Cu :   return true;
        case Fe :   return true;
        default :   return false;
    };
}
Material::MaterialType Material::Type() const {
    return type;
}
/****************************** class Body ******************************/
Body::BodyType Body::Type() const {
    return type;
}
Material Body::MaterialType() const {
    return material;
}
Body::Body(Body::BodyType t, Material m)
    : type(t)
    , material(m) {}
/**************************** class BoxBody ****************************/
BoxBody::BoxBody(double h, double l, double w, Material m)
	: Body(Body::Box, m)
{
    height = h;
    length = l;
    width = w;
}
double BoxBody::BasePerimeter() const {
    return 2*(length + width);
}
double BoxBody::BaseArea() const {
    return  length*width;
}
double BoxBody::SideArea() const {
    return BasePerimeter()*height;
}
double BoxBody::Height() const {
    return height;
}
double BoxBody::Length() const {
    return length;
}
double BoxBody::Width() const {
    return width;
}
double BoxBody::BaseIncircleRadius() const {
    return length == width ? length/2. : -1;
}
double BoxBody::BaseExcircleRadius() const {
    return sqrt(SQR(length)+SQR(width))/2.;
}
double BoxBody::Volume() const {
    return height*length*width;
}
double BoxBody::BasePerimeterStretchedBy(double eachSideAmount) const {
    return BasePerimeter() + 4*eachSideAmount;
}
/************************* class CylinderBody *************************/
CylindricBody::CylindricBody(double h, double r, Material m)
	: Body(Body::Cylinder, m)
{
    height = h;
    radius = r;
}
double CylindricBody::BasePerimeter() const {
    return 2*constants::pi*radius;
}
double CylindricBody::BaseArea() const {
    return constants::pi*SQR(radius);
}
double CylindricBody::SideArea() const {
    return BasePerimeter()*height;
}
double CylindricBody::Height() const {
    return height;
}
double CylindricBody::Radius() const {
    return radius;
}
double CylindricBody::BaseIncircleRadius() const {
    return radius;
}
double CylindricBody::BaseExcircleRadius() const {
    return radius;
}
double CylindricBody::Volume() const {
    return BaseArea()*height;
}
double CylindricBody::BasePerimeterStretchedBy(double radiusAmount) const {
    return BasePerimeter() + 2*constants::pi*radiusAmount;
}
/****************************** class Wire ******************************/
Wire::Wire(ElectricalEnvironment* e, double len, double metalD, double isolatedD, Material m)
{
    if (m.IsMetal()) {
        env = e;
        length = len;
        metalDiameter = metalD;
        isolatedDiameter = isolatedD;
	metal = m;
    }
}
double Wire::Length() const {
    return length;
}
double Wire::MetalDiameter() const {
    return metalDiameter;
}
double Wire::IsolatedDiameter() const {
    return isolatedDiameter;
}
double Wire::CutArea(bool isolated) const {
    return isolated ? constants::pi*SQR(isolatedDiameter/2.) : constants::pi*SQR(metalDiameter/2.);
}
double Wire::Mass(bool treatIsolationAsMetal) const {
    return CutArea(!treatIsolationAsMetal)*length*metal.Density();
}
double Wire::Resistance() const {
    return metal.Resistivity()*length/CutArea(false);
}
double Wire::Heating() const {
    return Resistance()*SQR(env->Current());
}
double Wire::TimeToMelt() const {
    double dm = CutArea()*IsolatedDiameter()*metal.Density();
    return metal.HeatingCapacity()*(metal.MeltingPoint() - env->Temperature())*dm/Heating();
}
Material Wire::Metal() const {
    return metal;
}
/****************************** class Coil ******************************/
Coil::Coil(ElectricalEnvironment* e, Wire w, PackingType pack, CylindricBody core) 
    : wire(w)
{
    env = e;
    packing = pack;
    coreBody = new CylindricBody(core);
    inductivity = 0.;
    winding();
}
Coil::Coil(ElectricalEnvironment* e, Wire w, PackingType pack, BoxBody core)
    : wire(w)
{
    env = e;
    packing = pack;
    coreBody = new BoxBody(core);
    inductivity = 0.;
    winding();
}
Coil::~Coil() {
    delete coreBody;
}
size_t Coil::Turns() const {
    return params.turns;
}
size_t Coil::Layers() const {
    return params.layers;
}
double Coil::WindingThickness() const {
    return params.thickness;
}
Coil::PackingType Coil::Packing() const {
   return packing; 
}
Material Coil::CoreMaterial() const {
    return coreBody->MaterialType();
}
double Coil::Inductivity() const {
    return inductivity;
}
double Coil::Resistance() const {
    return wire.Resistance() + env->Frequancy()*Inductivity();
}
double Coil::stretchAmount(double layer_iter) const {
   switch(packing) {
       case Quadratic :    return (layer_iter-1)*wire.IsolatedDiameter();
       case Compact   :    return (layer_iter-1)*wire.IsolatedDiameter()*constants::sqrt3/2.;
   };
}
void Coil::winding() {
    double full_layer_turns = coreBody->Height()/wire.IsolatedDiameter();
    double winded = 0.;
    double turns = 0., turn_perimeter = 0.;  
    unsigned int layer_iter = 0.;
    unsigned int iters = 0;
    while ( wire.Length() > winded ) {
        ++layer_iter;
	++iters;
        turn_perimeter = coreBody->BasePerimeterStretchedBy(stretchAmount(layer_iter));
        turns = (wire.Length() - winded)/turn_perimeter;
        winded += turn_perimeter*full_layer_turns;
        /** Inductivity **/
        double curr_turn_radius = 0., turn_radius = coreBody->BaseExcircleRadius() + stretchAmount(layer_iter);     // TODO: Do smth with rectangles
        double curr_distance = 0.;
	std::cout << "Turns on " << layer_iter << "layer - " << ( turns > full_layer_turns ? full_layer_turns : turns) << std::endl;
	std::cout << "\t\t\tIterations pass: " << iters << std::endl;
        for(size_t t = 1; t < ( turns > full_layer_turns ? full_layer_turns : turns); ++t) {
++iters;
            inductivity += (turns - t)*turnInductivity(t*wire.IsolatedDiameter(), turn_radius);                      // inductivity between current layers' turns
iters += 100;
            for(size_t l = layer_iter-1; l > 0; --l) {
++iters;
                curr_turn_radius = coreBody->BaseExcircleRadius() + stretchAmount(l);
                for(size_t prev_t = 0; prev_t < full_layer_turns; ++prev_t) {
++iters;
                    curr_distance = std::abs(t - prev_t)*wire.IsolatedDiameter();
                    inductivity += turnInductivity(curr_distance, curr_turn_radius, turn_radius);
iters+=100;
                }
            }
        }
        /**             **/
    }
    params.layers = layer_iter;
    params.turns = turns + full_layer_turns*(params.layers - 1);
    params.thickness = (params.layers - 1)*wire.IsolatedDiameter()*(Compact == packing ? constants::sqrt3/2. : 1) + wire.IsolatedDiameter();
}
double Coil::turnInductivity(double distanceBetweenTurns, double firstTurnRadius, double secondTurnRadius) {
    if( secondTurnRadius == -1.) {
        secondTurnRadius = firstTurnRadius;
    }
    double k = 2.*sqrt(firstTurnRadius*secondTurnRadius)/sqrt(SQR(firstTurnRadius + secondTurnRadius) + SQR(distanceBetweenTurns));
    aux::Function<double, double>* f1 = new aux::FunctionElipticIntegralOfFirstKind(k);
    aux::Function<double, double>* f2 = new aux::FunctionElipticIntegralOfSecondKind(k);

    double I1 = aux::CalculateDefineIntegralWithSimpsonMethod<double, double>(f1, 0., constants::pi/2., 100);
    double I2 = aux::CalculateDefineIntegralWithSimpsonMethod<double, double>(f2, 0., constants::pi/2., 100);
    
    delete f1;
    delete f2;

    return -constants::u0*sqrt(firstTurnRadius*secondTurnRadius)*((k - 2./k)*I1 + 2./k*I2);
}
}

using namespace CoilCalculation;

void __test() {
	ElectricalEnvironment env;
	CylindricBody core(.1, .05, Material::Fe);
	Wire wire(&env, 5164.97, .0005, .0005, Material::Cu);
	Coil coil(&env, wire, Coil::Quadratic, core);

	std::cout << "Winding turns    : " << coil.Turns() << std::endl;
	std::cout << "Winding layers   : " << coil.Layers() << std::endl;
	std::cout << "Winding thickness: " << coil.WindingThickness() << std::endl;
	std::cout << "!!! Inductivity  : " << coil.Inductivity() << std::endl;
	std::cout << "Coil Resistance  : " << coil.Resistance() << std::endl; 
}


int main() {
	__test();
	return 0;
}
