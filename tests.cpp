#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

// Functions to be tested are located here
#include "calculate.hpp"


TEST_CASE( "Calculate LJ", "[calculate]" ) {
    REQUIRE( calculate_LJ(1.0) == Approx(0.0) );
}

TEST_CASE( "Calculate Distance", "[calculate]" ) {
    REQUIRE( calculate_distance({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}) == Approx(0.0) );
    REQUIRE( calculate_distance({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}) == Approx(1.0) );
    REQUIRE( calculate_distance({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}) == Approx(1.0) );
    REQUIRE( calculate_distance({0.0, 1.0, 2.0}, {3.0, 4.0, 5.0}) == Approx(5.19615) );
}
