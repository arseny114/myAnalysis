# gaudi_add_header_only_library(ParticleIDLib)

# Modules
gaudi_add_module(myAnalysis
                 SOURCES src/myAnalysis.cpp
                 LINK GearSvc
                      SimplePIDSvc
                      Gaudi::GaudiAlgLib
                      Gaudi::GaudiKernel
                      k4FWCore::k4FWCore
                      DetInterface
                      DetSegmentation
                      DataHelperLib
                      ${GEAR_LIBRARIES} 
                      ${GSL_LIBRARIES} 
                      ${ROOT_LIBRARIES}
                      ${FASTJET_LIBRARY}
		                  EDM4CEPC::edm4cepc EDM4CEPC::edm4cepcDict
                      EDM4HEP::edm4hep EDM4HEP::edm4hepDict k4FWCore::k4FWCore
)

target_include_directories(myAnalysis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>/include
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(TARGETS myAnalysis
  EXPORT CEPCSWTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  COMPONENT dev)
