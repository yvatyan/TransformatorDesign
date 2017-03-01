#include <iostream>
#include <string>
#include <cmath>
#include <cstdio>
using namespace std;

#define SQR(x) ((x)*(x))

//TODO: Design following classes: BoxBody, CylindricBody, Wire, Winding, Coil

namespace constant {
    const double pi    = M_PI;
    const double u0    = 4*pi*1e-7;
    const double sqrt3 = sqrt(3.);
}

enum material {
    Cu,
    Fe
};
enum packingType {
    quadratic,
    compact
};

class Box {
    public:
        double height;
        double length;
        double width;
        Box(double h, double l, double w) {
            height = h;
            length = l;
            width  = w;
        }
        double perimeter() const {
            return 4*(height + length + width);
        }
        double widthSizePerimeter() const {
            return 2*(height + width);
        }
        double lengthSizePerimeter() const {
            return 2*(height + length);
        }
        double basePerimeter() const {
            return 2*(length + width);
        }
        double area() const {
            2*(height*length + length*width + width*height);
        }
        double widthSideArea() const {
            return height*width;
        }
        double lengthSideArea() const {
            return height*length;
        }
        double baseArea() const {
            return length*width;
        }
        double volume() const {
            return height*length*width;
        }
};

double metal_density (material metal) { // -> kg/m^3
    switch(metal) {
        case Cu :   return 8930.;
        case Fe :   return 7850.;
        default :   return -1;
    };
}

double metal_resistivity (material metal) { // -> Ohm*m [at 20C]
    switch(metal) {
        case Cu :   return 1.68e-8;
        case Fe :   return 9.71e-8;
        default :   return -1;
    };
}

double material_heating_capacity (material mat) { // -> J/(kg*K) 
     switch(mat) {
        case Cu :   return 385;
        case Fe :   return 444;
        default :   return -1;
    };
}

double material_permeability(material mat) {
    switch(mat) {
        case Cu :   return 1.26e-6;
        case Fe :   return 6.30e-3;
        default :   return -1;
    };
}

double wire_mass (double length, double diameter, material metal) { // m, m -> kg //TODO: Add to wire class
    double V = SQR(diameter/2.)*constant::pi*length;
    return V*metal_density(metal);
}

double heating (double mass, double temperature1, double temperature2, material mat) { // kg -> J //TODO: Merge with wire_heating
    return material_heating_capacity(mat)*(temperature2 - temperature1)*mass; 
}

double wire_resistance (double length, double diameter, material metal) { // m, m, -> Ohm   //TODO: Add to wire class
    double S = SQR(diameter/2.)*constant::pi;
    return metal_resistivity(metal)*length/S; 
}

double winding (double wireLength, double wireDiameter, Box coreSize, packingType packing = quadratic) { // m, m, (m, m, m) -> turns //TODO: Make it object Winding{turns, packing, wire{length, diameter, isoDiameter, material}, bodyShape{cylindric{radius, height}|box{height, width, lebgth}}, coreMaterial}
    double k = coreSize.height/wireDiameter; 
    double Sum = 0.;
    double N_layer = 0.;  
    unsigned int iter = 0.;
    while ( wireLength > Sum ) {
    std::cout << "N_layer: " << N_layer << std::endl;
        ++iter;
        if ( quadratic == packing ) {
            N_layer = (wireLength - Sum)/(coreSize.basePerimeter() + 4*iter*wireDiameter);
            Sum += (coreSize.basePerimeter() + 4*iter*wireDiameter)*k;
        } else if ( compact == packing ) {
            N_layer = (wireLength - Sum)/(coreSize.basePerimeter() + 4*wireDiameter + 4*(iter-1)*constant::sqrt3*wireDiameter/2.);
            Sum += (coreSize.basePerimeter() + 4*wireDiameter + 4*(iter-1)*constant::sqrt3*wireDiameter/2.)*k;
        }
        std::cout << "Sum: " << Sum << std::endl;
    }
    return N_layer + k*(iter-1);
}

double coil_inductivity (double wireLength, double wireDiameter, Box coreSize, material coreMetal) { // m, m, (m, m, m) -> H    //TODO: Add to coil class
    double N = winding(wireLength, wireDiameter, coreSize);
    double W = ((int)(N/(coreSize.height/wireDiameter)) + 1)*wireDiameter; // TODO: This is for quadratic packing, Make it for compact class
    double R = W/2 + sqrt(coreSize.baseArea() / constant::pi);
    std::cout << "Packing Thichness: " << W << std::endl;
    return (/*material_permeability(coreMetal)*/constant::u0*SQR(N)*SQR(R))/(2*constant::pi*(6*R + 9*coreSize.height + 10*W)); 
}

double coil_resistance (double wireLength, double wireDiameter, Box coilSize, double frequency, material coreMetal, material wireMetal) { // m, m, (m, m, m), Hz -> Ohm //TODO: Add to coil class
    return wire_resistance(wireLength, wireDiameter, wireMetal) + frequency*coil_inductivity(wireLength, wireDiameter, coilSize, coreMetal);
}

double wire_heating (double resistance, double current) { // Ohm, A -> J/sec    // TODO: Add to wire class
    return resistance*SQR(current);
}

void _test() {
    double wireLength = 1.44;
    double wireDiameter = .001;
    Box coilSize(.01, .01, .01);
    material mat = Cu;
    std::cout << "Winding turns (quadratic): " << winding(wireLength, wireDiameter, coilSize) << std::endl;
    std::cout << "Winding turns (compact)  : " << winding(wireLength, wireDiameter, coilSize, compact) << std::endl;
    std::cout << "Inductivity: "   << 1000000*coil_inductivity(wireLength, wireDiameter, coilSize, mat) << std::endl;
    std::cout << "Resistance: "   << coil_resistance(wireLength, wireDiameter, coilSize, 0, Fe, Cu) << std::endl;
}

int main() {
    _test();
    return 0;
}
