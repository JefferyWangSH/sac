#ifndef SAC_INITIALIZER_H
#define SAC_INITIALIZER_H
#pragma once

/**
  *  This header file defines the interface class `SAC::Initializer`
  *  for the initialization of the relevant modules used in the SAC simulation.
  *  Program options are read from the toml configuration file
  *  using the static member function SAC::Initializer::toml_parse().
  */



namespace SAC {

    // ------------------------------------  SAC::Initializer class  --------------------------------------
    class Initializer {

        public:

            static void toml_parse();

            static void initial();



    };

} // namespace SAC

#endif // SAC_INITIALIZER_H
