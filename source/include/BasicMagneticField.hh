

#ifndef BasicMagneticField_H
#define BasicMagneticField_H

#include "G4UniformMagField.hh"

class G4FieldManager;



class BasicMagneticField: public G4UniformMagField
{
  public:
  
   BasicMagneticField(G4ThreeVector);  //  The value of the field
   BasicMagneticField();               //  A zero field
  ~BasicMagneticField();  
      
   //Set the field (0,0,fieldValue)
   void SetFieldValue(G4double fieldValue);
   void SetFieldValue(G4ThreeVector fieldVector);

protected:

  // Find the global Field Manager
  G4FieldManager* GetGlobalFieldManager();   // static 
};

#endif
