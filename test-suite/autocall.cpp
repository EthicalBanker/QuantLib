/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014-6 Catley Lakeman Ltd. 
 The moral rights of the author of the original version of this file, David Rees, have been asserted.

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "autocall.hpp"
#include "utilities.hpp"
#include <ql/quantlib.hpp>

#include <ql/experimental/exoticoptions/mcautocallengine.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/pricingengines/barrier/analyticbarrierengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/time/daycounters/actual360.hpp>

using namespace boost::unit_test_framework;

namespace {
    const static double cSamplePower = 10.0;

    void autocall_boost_check( 
        const std::string& method, 
        const double autocall_NPV, 
        const double expected,
        const double tolerance )
    {
        const double error = std::fabs( autocall_NPV - expected );

        if( error > tolerance)
            BOOST_ERROR(std::setprecision(16) << "\n" << method << " " \
                        << "autocall NPV: " << autocall_NPV << "\n"
                        << "expected:       = " << expected << "\n"
                        << "error:        = " << error );
    }
} // namespace

test_suite* AutocallTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Autocall tests");

    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondBeforeIssueBeforeStrike) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondBeforeIssueAfterStrike) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondOnCouponDate) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondOffCouponDate) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondAtMaturity) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBondAfterMaturity) );

    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBeforeIssueBeforeStrike) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testBeforeIssueAfterStrike) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testOnCouponDate) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testOffCouponDate) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testAtMaturity) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testAfterMaturity) );
    
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testExample) );
    
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsZeroStrikeCall) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsCallOption) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsDigitalOption) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsCallSpread) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsBarrierOption) );
    suite->add( QUANTLIB_TEST_CASE(&AutocallTest::testVsEverest) );

    return suite;
}

void AutocallTest::testBondBeforeIssueBeforeStrike() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date before issue and before strike" );   
    testBondStyle( "Testing BondStyle BeforeIssueBeforeStrike", BeforeIssueBeforeStrike, 1.0e-6 );
}

void AutocallTest::testBondBeforeIssueAfterStrike() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date before issue and after strike" );   
    testBondStyle( "Testing BondStyle BeforeIssueAfterStrike", BeforeIssueAfterStrike, 1.0e-6 );
}

void AutocallTest::testBondOnCouponDate() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date on coupon date" );   
    testBondStyle( "Testing BondStyle OnCouponDate", OnCouponDate, 1.0e-6 );
}

void AutocallTest::testBondOffCouponDate() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date not on a coupon date" );   
    testBondStyle( "Testing BondStyle OffCouponDate", OffCouponDate, 1.0e-6 );
}

void AutocallTest::testBondAtMaturity() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date on maturity" );   
    testBondStyle( "Testing BondStyle AtMaturity", AtMaturity, 1.0e-6 );
}

void AutocallTest::testBondAfterMaturity() {
    BOOST_TEST_MESSAGE("Testing Autocall option vs bond with " <<
                        "evaluation date after maturity" );   
    testBondStyle( "Testing BondStyle AfterMaturity", AfterMaturity, 1.0e-6 );
}

void AutocallTest::testBeforeIssueBeforeStrike() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date before issue and before strike" );   
    const Real expectedResult = 258049.12493518894;
    testNormal( "Testing BeforeIssueBeforeStrike", BeforeIssueBeforeStrike, expectedResult, 1.0e-6 );
}

void AutocallTest::testBeforeIssueAfterStrike() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date before issue and after strike" );   
    const Real expectedResult = 252116.10118340506;
    testNormal( "Testing BeforeIssueAfterStrike", BeforeIssueAfterStrike, expectedResult, 1.0e-6 );
}

void AutocallTest::testOnCouponDate() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date on coupon date" );   
    const Real expectedResult = 357751.18696049560;
    testNormal( "Testing OnCouponDate", OnCouponDate, expectedResult, 1.0e-6 );
}

void AutocallTest::testOffCouponDate() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date not on a coupon date" );   
    const Real expectedResult = 361908.97569314751;
    testNormal( "Testing OffCouponDate", OffCouponDate, expectedResult, 1.0e-6 );
}

void AutocallTest::testAtMaturity() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date on maturity" );   
    const Real expectedResult = 0.0;
    testNormal( "Testing AtMaturity", AtMaturity, expectedResult, 1.0e-6 );
}

void AutocallTest::testAfterMaturity() {
    BOOST_TEST_MESSAGE("Testing Autocall option with " <<
                        "evaluation date after maturity" );   
    const Real expectedResult = 0.0;
    testNormal( "Testing AfterMaturity", AfterMaturity, expectedResult, 1.0e-6 );
}

// A complete example of creating and calling npv.
void AutocallTest::testExample() {
    BOOST_TEST_MESSAGE("Testing Autocall option example" );
    
    // Create an autocall
    const Date today( 28, February, 2014 );
    const Date issueDate = today + 30;

    const Integer yearsToMaturity = 6;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1,1.0e6);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	const std::vector<Real> fixedCoupons( yearsToMaturity, 0.03 );
    const Schedule& callDates = fixedSchedule;
    std::vector<Real> callLevels( yearsToMaturity, 1.0 );
    for( Size i=0 ; i<yearsToMaturity ; ++i ) 
        callLevels[i] = 1.0 + i/100.0;    
    const std::vector<Real> callPayments( yearsToMaturity, 0.05 );
    
    const Date strikeDate = issueDate - 10;
    std::vector<Real> strikeLevels;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  QL_MAX_REAL;
    const Real redemptionUpsideGearing   =  1.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  1.0;
    const Real redemptionBarrier         =  0.6;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier 
    );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo 
    );
    
    // Create a dummy yield curve for discounting
    DayCounter dc = Actual360();
    Handle<YieldTermStructure> riskFreeRate(flatRate(today, 0.05, dc));

    // Create stochastic processes for the underlying indices
    std::vector<boost::shared_ptr<StochasticProcess1D> > processes(4);
    processes[0] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(
              Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(100.0))),
              Handle<YieldTermStructure>(flatRate(today, 0.01, dc)),
              riskFreeRate,
              Handle<BlackVolTermStructure>(flatVol(today, 0.30, dc))));
    processes[1] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(
              Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(110.0))),
              Handle<YieldTermStructure>(flatRate(today, 0.05, dc)),
              riskFreeRate,
              Handle<BlackVolTermStructure>(flatVol(today, 0.35, dc))));
    processes[2] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(
              Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(90.0))),
              Handle<YieldTermStructure>(flatRate(today, 0.04, dc)),
              riskFreeRate,
              Handle<BlackVolTermStructure>(flatVol(today, 0.25, dc))));
    processes[3] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(
              Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(105.0))),
              Handle<YieldTermStructure>(flatRate(today, 0.03, dc)),
              riskFreeRate,
              Handle<BlackVolTermStructure>(flatVol(today, 0.20, dc))));

    // Create a dummy correlation matrix
    Matrix correlation(4,4);
    correlation[0][0] = 1.00;
                    correlation[0][1] = 0.50;
                                    correlation[0][2] = 0.30;
                                                    correlation[0][3] = 0.10;
    correlation[1][0] = 0.50;
                    correlation[1][1] = 1.00;
                                    correlation[1][2] = 0.20;
                                                    correlation[1][3] = 0.40;
    correlation[2][0] = 0.30;
                    correlation[2][1] = 0.20;
                                    correlation[2][2] = 1.00;
                                                    correlation[2][3] = 0.60;
    correlation[3][0] = 0.10;
                    correlation[3][1] = 0.40;
                                    correlation[3][2] = 0.60;
                                                    correlation[3][3] = 1.00;
    
    // Create a process array
    boost::shared_ptr<StochasticProcessArray> process_array(
                          new StochasticProcessArray(processes, correlation));

    // Set Monte Carlo inputs
    const Size fixedSamples = 1023;
    const BigNatural seed = 86421;
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    // Calculate the net present value, and compare to a stored value
    const Real autocall_NPV = autocall.NPV();
    const Real expectedResult = 391693.9037542623;
    const Real tolerance = 1.0e-6;

    autocall_boost_check( "Test Example", autocall_NPV, expectedResult, tolerance );
}

// Test a bond style autocall 
void AutocallTest::testBondStyle(
    const std::string& method,
    const EvaluationDateType& evaluationDateType,
    const Real& tolerance )
{
    const Date todayDate( 4, July, 2014 );
    const Date issueDate = todayDate;

    boost::shared_ptr<Autocall> bondStyleAutocall;
    boost::shared_ptr<FixedRateBond> bond;
    const Size numUnderlyings = 1;
    const bool strikeLevelsSet = evaluationDateType != BeforeIssueBeforeStrike;

    // Create the autocall and equivalent bond
    createBondStyleAutocall( 
        issueDate,
        numUnderlyings,
        strikeLevelsSet,
        bondStyleAutocall, 
        bond );
  
    // Set the valuation date
    Date evaluationDate;
    switch( evaluationDateType ) {
    case BeforeIssueBeforeStrike :
        evaluationDate = todayDate - 30;
        break;
    case BeforeIssueAfterStrike :
        evaluationDate = todayDate - 5;
        break;
    case OnCouponDate :
        evaluationDate = bond->cashflows()[4]->date();
        break;
    case OffCouponDate :
        evaluationDate = bond->cashflows()[4]->date() + 45;
        break;
    case AtMaturity:
        evaluationDate = bond->maturityDate();
        break;
    case AfterMaturity :
        evaluationDate = bond->maturityDate() + 100;
        break;
    default :
        BOOST_FAIL( "Unknown valuationDateType" );
        break;
    }

    Date originalEvaluationDate = Settings::instance().evaluationDate(); 
    Settings::instance().evaluationDate() = evaluationDate;

    // Create market data
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketData(
        evaluationDate,
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = 1023;

    bondStyleAutocall->setPricingEngine(MakeMCAutocallEngine<PseudoRandom>(process_array)
                        .withStepsPerYear(10)
						.withSamples(fixedSamples)
						.withSeed(seed));

    boost::shared_ptr<PricingEngine> bondEngine(new DiscountingBondEngine(riskFreeRate));
    bond->setPricingEngine( bondEngine );

    const Real autocall_NPV = bondStyleAutocall->NPV();
    const Real expectedResult = bond->NPV();
    
    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// Test a normal autocall 
void AutocallTest::testNormal(
    const std::string& method,
    const EvaluationDateType& evaluationDateType,
    const Real& expectedResult,
    const Real& tolerance )
{
    const Date todayDate( 17, September, 2013 );
    const Date issueDate = todayDate;

    boost::shared_ptr<Autocall> normalAutocall;
    const Size numUnderlyings = 3;
    const bool strikeLevelsSet = evaluationDateType != BeforeIssueBeforeStrike;

    // Create the autocall and equivalent bond
    createNormalAutocall( 
        issueDate,
        numUnderlyings,
        strikeLevelsSet,
        normalAutocall );
  
    // Set the valuation date
    Date evaluationDate;
    switch( evaluationDateType ) {
    case BeforeIssueBeforeStrike :
        evaluationDate = todayDate - 30;
        break;
    case BeforeIssueAfterStrike :
        evaluationDate = todayDate - 5;
        break;
    case OnCouponDate :
        evaluationDate = normalAutocall->getFixedSchedule()[3];
        break;
    case OffCouponDate :
        evaluationDate = normalAutocall->getFixedSchedule()[3] + 45;
        break;
    case AtMaturity:
        evaluationDate = normalAutocall->getFixedSchedule().endDate();
        break;
    case AfterMaturity :
        evaluationDate = normalAutocall->getFixedSchedule().endDate() + 100;
        break;
    default :
        BOOST_FAIL( "Unknown valuationDateType" );
        break;
    }

    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = evaluationDate;

    // Create market data
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketData(
        evaluationDate,
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = 1023;

    normalAutocall->setPricingEngine(MakeMCAutocallEngine<PseudoRandom>(process_array)
                        .withStepsPerYear(10)
						.withSamples(fixedSamples)
						.withSeed(seed));

    const Real autocall_NPV = normalAutocall->NPV();
    
    Settings::instance().evaluationDate() = originalEvaluationDate;    
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// Create an autocall with same cashflows as fixed bond
void AutocallTest::createBondStyleAutocall(
    const Date& issueDate,
    const Size& numUnderlyings,
    const bool strikeLevelsSet,
    boost::shared_ptr<Autocall>& bondStyleAutocall,
    boost::shared_ptr<FixedRateBond>& bond
    ) 
{
    const Integer yearsToMaturity = 6;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1,1.0e6);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	const std::vector<Real> fixedCoupons( yearsToMaturity, 0.03 );
    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.05 );

    const Date strikeDate = issueDate - 10;
    std::vector<Real> strikeLevels;
    if( strikeLevelsSet ) {
        for( Size i=0 ; i<numUnderlyings ; ++i )
            strikeLevels.push_back( 100.0 + i );
    }

    const Real redemptionCash            = 1.0;
    const Real redemptionCap             = 0.0;
    const Real redemptionUpsideGearing   = 0.0;
    const Real redemptionFloor           = 0.0;
    const Real redemptionCap2            = 0.0;
    const Real redemptionDownsideGearing = 0.0;
    const Real redemptionBarrier         = QL_MAX_REAL;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
    );

    bondStyleAutocall = boost::shared_ptr<Autocall>( new Autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    ) );

    // Create an equivalent bond
    const Natural settlementDays = 0;
    const DayCounter dc = SimpleDayCounter();
    bond = boost::shared_ptr<FixedRateBond>( new FixedRateBond(
        settlementDays,
        notionals[0],
        fixedSchedule,
        fixedCoupons,
        dc 
    ) );
}

// Create a normal autocall
void AutocallTest::createNormalAutocall(
    const Date& issueDate,
    const Size& numUnderlyings,
    const bool strikeLevelsSet,
    boost::shared_ptr<Autocall>& autocall
    ) 
{
    const Integer yearsToMaturity = 6;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1,1.0e6);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	const std::vector<Real> fixedCoupons( yearsToMaturity, 0.03 );
    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, 1.0 );
    const std::vector<Real> callPayments( yearsToMaturity, 0.05 );

    const Date strikeDate = issueDate - 10;
    std::vector<Real> strikeLevels;
    if( strikeLevelsSet ) {
        for( Size i=0 ; i<numUnderlyings ; ++i )
            strikeLevels.push_back( 100.0 + i );
    }

    const Real redemptionCash            = 1.0;
    const Real redemptionCap             = 2.0;
    const Real redemptionUpsideGearing   = 0.9;
    const Real redemptionFloor           = 0.5;
    const Real redemptionCap2            = 3.0;
    const Real redemptionDownsideGearing = 1.1;
    const Real redemptionBarrier         = 0.6;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    autocall = boost::shared_ptr<Autocall>( new Autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    ) );
}

// Create dummy market data as of valuation date
void AutocallTest::createMarketData(
    const Date& valuationDate,
    const Size& numUnderlyings,
    Handle<YieldTermStructure>& riskFreeRate,
    boost::shared_ptr<StochasticProcessArray>& process_array
    )
{
    // Create a flat yield curve
    const DayCounter dc = Actual365Fixed();
    riskFreeRate = Handle<YieldTermStructure>( flatRate( valuationDate, 0.05, dc ) );

    // Create one process per underlying
    std::vector<boost::shared_ptr<StochasticProcess1D> > processes( numUnderlyings );
    for( Size i=0 ; i<numUnderlyings ; ++i ) {
        processes[i] = boost::shared_ptr<StochasticProcess1D>(
			new BlackScholesMertonProcess(
				  Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(100.0+i))),
				  Handle<YieldTermStructure>(flatRate(valuationDate, i/100.0, dc)),
				  riskFreeRate,
				  Handle<BlackVolTermStructure>(flatVol(valuationDate, 0.05+i/20.0, dc))));
    }

    // Create correlation matrix
    Matrix correlation( numUnderlyings, numUnderlyings );
    for( Size i=0 ; i<numUnderlyings ; ++i ) {
        for( Size j=0 ; j<numUnderlyings ; ++j ) {
            Real ratio = static_cast<double>(i+1)/static_cast<double>(j+1);
            correlation[i][j] = i>j ? 1.0/ratio : ratio;
        }
    }

    process_array = boost::shared_ptr<StochasticProcessArray>(
							  new StochasticProcessArray(processes, correlation));

}

// Create some more realistic market data as of valuation date
void AutocallTest::createMarketDataRealistic(
    const Size& numUnderlyings,
    Handle<YieldTermStructure>& riskFreeRate,
    boost::shared_ptr<StochasticProcessArray>& process_array
    )
{
    const Date today(1, July, 2014);
    
    // Create a yield curve from actual market data
    const DayCounter dc = Actual365Fixed();
    
    std::vector<Date> dates;
    dates.push_back( today             );
    dates.push_back( Date( 3,Jul,2014) );
    dates.push_back( Date( 7,Jul,2014) );
    dates.push_back( Date(10,Jul,2014) );
    dates.push_back( Date(17,Jul,2014) );
    dates.push_back( Date( 1,Aug,2014) );
    dates.push_back( Date( 1,Sep,2014) );
    dates.push_back( Date( 1,Oct,2014) );
    dates.push_back( Date( 1,Nov,2014) );
    dates.push_back( Date( 3,Dec,2014) );
    dates.push_back( Date( 1,Jan,2015) );
    dates.push_back( Date( 3,Feb,2015) );
    dates.push_back( Date( 3,Mar,2015) );
    dates.push_back( Date( 1,Apr,2015) );
    dates.push_back( Date( 4,May,2015) );
    dates.push_back( Date( 3,Jun,2015) );
    dates.push_back( Date( 1,Jul,2015) );
    dates.push_back( Date( 5,Oct,2015) );
    dates.push_back( Date( 4,Jan,2016) );
    dates.push_back( Date( 4,Apr,2016) );
    dates.push_back( Date( 1,Jul,2016) );
    dates.push_back( Date( 3,Oct,2016) );
    dates.push_back( Date( 3,Jan,2017) );
    dates.push_back( Date( 3,Apr,2017) );
    dates.push_back( Date( 1,Jul,2017) );
    dates.push_back( Date( 3,Oct,2017) );
    dates.push_back( Date( 3,Jan,2018) );
    dates.push_back( Date( 3,Apr,2018) );
    dates.push_back( Date( 1,Jul,2018) );
    dates.push_back( Date( 3,Oct,2018) );
    dates.push_back( Date( 3,Jan,2019) );
    dates.push_back( Date( 3,Apr,2019) );
    dates.push_back( Date( 1,Jul,2019) );
    dates.push_back( Date( 3,Oct,2019) );
    dates.push_back( Date( 3,Jan,2020) );
    dates.push_back( Date( 3,Apr,2020) );
    dates.push_back( Date( 1,Jul,2020) );
    dates.push_back( Date( 5,Oct,2020) );
    dates.push_back( Date( 4,Jan,2021) );
    dates.push_back( Date( 5,Apr,2021) );
    dates.push_back( Date( 5,Jul,2021) );
    dates.push_back( Date( 4,Oct,2021) );
    dates.push_back( Date( 3,Jan,2022) );
    dates.push_back( Date( 4,Apr,2022) );
    dates.push_back( Date( 4,Jul,2022) );
    dates.push_back( Date( 3,Oct,2022) ); 
    dates.push_back( Date( 3,Jan,2023) ); 
    dates.push_back( Date( 3,Apr,2023) ); 
    dates.push_back( Date( 3,Jul,2023) ); 
    dates.push_back( Date( 3,Oct,2023) ); 
    dates.push_back( Date( 3,Jan,2024) ); 
    dates.push_back( Date( 3,Apr,2024) ); 
    dates.push_back( Date( 3,Jul,2024) ); 
    dates.push_back( Date( 3,Jul,2025) ); 
    dates.push_back( Date( 3,Jul,2026) ); 
    dates.push_back( Date( 5,Jul,2027) ); 
    dates.push_back( Date( 3,Jul,2028) ); 
    dates.push_back( Date( 3,Jul,2029) ); 
    dates.push_back( Date( 3,Jul,2030) ); 
    dates.push_back( Date( 3,Jul,2031) ); 
    dates.push_back( Date( 5,Jul,2032) ); 
    dates.push_back( Date( 4,Jul,2033) ); 
    dates.push_back( Date( 3,Jul,2034) );
     
    std::vector<DiscountFactor> discountFactor;
    discountFactor.push_back( 1.0      ); 
    discountFactor.push_back( 0.999974 );
    discountFactor.push_back( 0.999921 );
    discountFactor.push_back( 0.999881 );
    discountFactor.push_back( 0.999787 );
    discountFactor.push_back( 0.999577 );
    discountFactor.push_back( 0.999105 );
    discountFactor.push_back( 0.998593 );
    discountFactor.push_back( 0.997941 );
    discountFactor.push_back( 0.997185 );
    discountFactor.push_back( 0.996421 );
    discountFactor.push_back( 0.995539 );
    discountFactor.push_back( 0.994779 );
    discountFactor.push_back( 0.993821 );
    discountFactor.push_back( 0.992818 );
    discountFactor.push_back( 0.991865 );
    discountFactor.push_back( 0.990859 );
    discountFactor.push_back( 0.987179 );
    discountFactor.push_back( 0.982984 );
    discountFactor.push_back( 0.978283 );
    discountFactor.push_back( 0.973249 );
    discountFactor.push_back( 0.967709 );
    discountFactor.push_back( 0.961862 );
    discountFactor.push_back( 0.955661 );
    discountFactor.push_back( 0.949191 );
    discountFactor.push_back( 0.942662 );
    discountFactor.push_back( 0.935976 );
    discountFactor.push_back( 0.929039 );
    discountFactor.push_back( 0.922035 );
    discountFactor.push_back( 0.915039 );
    discountFactor.push_back( 0.908029 );
    discountFactor.push_back( 0.900891 );
    discountFactor.push_back( 0.893736 );
    discountFactor.push_back( 0.886499 );
    discountFactor.push_back( 0.879317 );
    discountFactor.push_back( 0.872038 );
    discountFactor.push_back( 0.864832 );
    discountFactor.push_back( 0.857254 );
    discountFactor.push_back( 0.849897 );
    discountFactor.push_back( 0.842384 );
    discountFactor.push_back( 0.835000 );
    discountFactor.push_back( 0.827643 );
    discountFactor.push_back( 0.820300 );
    discountFactor.push_back( 0.812829 );
    discountFactor.push_back( 0.805373 );
    discountFactor.push_back( 0.798061 );
    discountFactor.push_back( 0.790662 );
    discountFactor.push_back( 0.783271 );
    discountFactor.push_back( 0.775859 );
    discountFactor.push_back( 0.768570 );
    discountFactor.push_back( 0.761273 );
    discountFactor.push_back( 0.753979 );
    discountFactor.push_back( 0.746705 );
    discountFactor.push_back( 0.718427 );
    discountFactor.push_back( 0.690457 );
    discountFactor.push_back( 0.663799 );
    discountFactor.push_back( 0.637914 );
    discountFactor.push_back( 0.612633 );
    discountFactor.push_back( 0.589227 );
    discountFactor.push_back( 0.566584 );
    discountFactor.push_back( 0.544622 );
    discountFactor.push_back( 0.523604 );
    discountFactor.push_back( 0.503429 );
    
    riskFreeRate = Handle<YieldTermStructure>( 
        boost::shared_ptr<YieldTermStructure>( 
            new InterpolatedDiscountCurve<LogLinear>(dates,discountFactor,dc) ) );

    // Create one process per underlying
    std::vector<boost::shared_ptr<StochasticProcess1D> > processes( numUnderlyings );

    // Add in some realistic FTSE data, for underlying 0
    std::vector<Real> strikes;
    strikes.push_back( 5442.34 );
    strikes.push_back( 5782.48 ); 
    strikes.push_back( 6122.63 ); 
    strikes.push_back( 6462.77 );
    strikes.push_back( 6802.92 );
    strikes.push_back( 7143.07 ); 
    strikes.push_back( 7483.21 ); 
    strikes.push_back( 7823.36 ); 
    strikes.push_back( 8163.50 ); 
    strikes.push_back( 8503.65 );

    std::vector<Date> expirations;
    expirations.push_back( Date(01, Jul, 2015) );
    expirations.push_back( Date(01, Jul, 2016) ); 
    expirations.push_back( Date(01, Jul, 2017) );
    expirations.push_back( Date(01, Jul, 2018) ); 
    expirations.push_back( Date(01, Jul, 2019) );
    expirations.push_back( Date(01, Jul, 2020) );
    
    Matrix volMatrix( strikes.size(), expirations.size() );
    
    //5442.34
    volMatrix[0][0] = .1796;
    volMatrix[0][1] = .1848;
    volMatrix[0][2] = .1909;
    volMatrix[0][3] = .1957;
    volMatrix[0][4] = .1990;
    volMatrix[0][5] = .2017;
    
    //5782.48
    volMatrix[1][0] = .1646;
    volMatrix[1][1] = .1738;
    volMatrix[1][2] = .1819;
    volMatrix[1][3] = .1879;
    volMatrix[1][4] = .1921;
    volMatrix[1][5] = .1953;
    
    //6122.63
    volMatrix[2][0] = .1503;
    volMatrix[2][1] = .1636;
    volMatrix[2][2] = .1736;
    volMatrix[2][3] = .1807;
    volMatrix[2][4] = .1857;
    volMatrix[2][5] = .1895;
    
    //6462.77
    volMatrix[3][0] = .1368;
    volMatrix[3][1] = .1541;
    volMatrix[3][2] = .1659;
    volMatrix[3][3] = .1741;
    volMatrix[3][4] = .1798;
    volMatrix[3][5] = .1842;
    
    //6802.92
    volMatrix[4][0] = .1244;
    volMatrix[4][1] = .1455;
    volMatrix[4][2] = .1590;
    volMatrix[4][3] = .1682;
    volMatrix[4][4] = .1745;
    volMatrix[4][5] = .1793;
    
    //7143.07
    volMatrix[5][0] = .1139;
    volMatrix[5][1] = .1381;
    volMatrix[5][2] = .1529;
    volMatrix[5][3] = .1628;
    volMatrix[5][4] = .1697;
    volMatrix[5][5] = .1750;
    
    //7483.21
    volMatrix[6][0] = .1061;
    volMatrix[6][1] = .1319;
    volMatrix[6][2] = .1476;
    volMatrix[6][3] = .1582;
    volMatrix[6][4] = .1655;
    volMatrix[6][5] = .1710;
    
    //7823.36
    volMatrix[7][0] = .1017;
    volMatrix[7][1] = .1271;
    volMatrix[7][2] = .1432;
    volMatrix[7][3] = .1541;
    volMatrix[7][4] = .1618;
    volMatrix[7][5] = .1675;
    
    //8163.50
    volMatrix[8][0] = .1005;
    volMatrix[8][1] = .1238;
    volMatrix[8][2] = .1397;
    volMatrix[8][3] = .1508;
    volMatrix[8][4] = .1586;
    volMatrix[8][5] = .1645;
    
    //8503.65
    volMatrix[9][0] = .1013;
    volMatrix[9][1] = .1219;
    volMatrix[9][2] = .1372;
    volMatrix[9][3] = .1480;
    volMatrix[9][4] = .1559;
    volMatrix[9][5] = .1618;
    
    const Date evaluationDate = today;
    const Calendar calendar = UnitedKingdom(UnitedKingdom::Exchange);
    const DayCounter dayCounter = ActualActual();

    const boost::shared_ptr<BlackVarianceSurface> volatilitySurface =
        boost::make_shared<BlackVarianceSurface>(
                evaluationDate, calendar, expirations, strikes, volMatrix, dayCounter );
    volatilitySurface->setInterpolation<Bicubic>();
    volatilitySurface->enableExtrapolation(true);

    /* Realistic surface doesn't work for some reason */
    processes[0] = boost::shared_ptr<StochasticProcess1D>(
        new BlackScholesMertonProcess(
            Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(6802.92))),
            Handle<YieldTermStructure>(flatRate(evaluationDate, 0.035, dayCounter)),
            riskFreeRate,
            Handle<BlackVolTermStructure>(volatilitySurface)));
    

    for( Size i=0 ; i<numUnderlyings ; ++i ) {
        processes[i] = boost::shared_ptr<StochasticProcess1D>(
			new BlackScholesMertonProcess(
				  Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(6802.92+100+i))),
				  Handle<YieldTermStructure>(flatRate(evaluationDate, 0.035+i/100.0, dc)),
				  riskFreeRate,
				  Handle<BlackVolTermStructure>(flatVol(evaluationDate, 0.05+i/20.0, dc))));
    }
    
    // Create correlation matrix
    Matrix correlation( numUnderlyings, numUnderlyings );
    for( Size i=0 ; i<numUnderlyings ; ++i ) {
        for( Size j=0 ; j<numUnderlyings ; ++j ) {
            Real ratio = static_cast<double>(i+1)/static_cast<double>(j+1);
            correlation[i][j] = i>j ? 1.0/ratio : ratio;
        }
    }

    process_array = boost::shared_ptr<StochasticProcessArray>(
							  new StochasticProcessArray(processes, correlation));
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs a zero strike call option
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsZeroStrikeCall()
{
    static const std::string method = "Testing Autocall vs Zero Strike Call";
    BOOST_TEST_MESSAGE( method );

    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	const std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.0 );

    const Date strikeDate = issueDate;
    std::vector<Real> spotPrices(1, 6710.45);
    std::vector<Real> strikeLevels = spotPrices;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  QL_MAX_REAL;
    const Real redemptionUpsideGearing   =  1.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  1.0;
    const Real redemptionBarrier         =  QL_MAX_REAL;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create a european option with zero strake
    boost::shared_ptr<PlainVanillaPayoff> payoff(new
        PlainVanillaPayoff( Option::Call, 0.0 ) );

    const Date exDate = maturityDate;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));
    EuropeanOption zeroStrikeCall( payoff, exercise );

    // Create market data
    const Size numUnderlyings = 1;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketDataRealistic(
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    boost::shared_ptr<GeneralizedBlackScholesProcess> process =
        boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(
                                                      process_array->process(0));
    QL_REQUIRE(process, "Black-Scholes process required");

    /* Could test vs Monte Carl
    zeroStrikeCall.setPricingEngine(MakeMCEuropeanEngine<LowDiscrepancy>(process)
                        .withAntitheticVariate(true)
                        .withStepsPerYear(stepsPerYear)
                        .withSamples(fixedSamples)
                        .withSeed(seed));
    */
    zeroStrikeCall.setPricingEngine( 
        boost::shared_ptr<PricingEngine>( new AnalyticEuropeanEngine(process) ) 
        );

    const Real autocall_NPV = autocall.NPV();
    const Real option_NPV = zeroStrikeCall.NPV();
    const Real expectedResult = option_NPV/spotPrices[0];
    const Real tolerance = 1.0e-2;
   
    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs a call option
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsCallOption() {
    static const std::string method = "Testing Autocall vs Call Option";
    BOOST_TEST_MESSAGE( method );

    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    fixedCoupons.back() = -1.0;

    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.0 );

    const Date strikeDate = issueDate;
    std::vector<Real> spotPrices( 1, 6710.45 );
    std::vector<Real> strikeLevels = spotPrices;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  QL_MAX_REAL;
    const Real redemptionUpsideGearing   =  1.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  0.0;
    const Real redemptionBarrier         =  1.0;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create a european at the money option
    boost::shared_ptr<PlainVanillaPayoff> payoff(new
            PlainVanillaPayoff( Option::Call, spotPrices[0] ) );

    const Date exDate = maturityDate;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));
    EuropeanOption callOption( payoff, exercise );

    // Create market data
    const Size numUnderlyings = 1;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketDataRealistic(
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const BigNatural stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    boost::shared_ptr<GeneralizedBlackScholesProcess> process =
        boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(
                                                      process_array->process(0));
    QL_REQUIRE(process, "Black-Scholes process required");

    callOption.setPricingEngine(
        boost::shared_ptr<PricingEngine>( new AnalyticEuropeanEngine(process) ) 
        );

    const Real autocall_NPV = autocall.NPV();
    const Real option_NPV = callOption.NPV();
    const Real expectedResult = option_NPV/spotPrices[0];
    const Real tolerance = 1.0e-2;

    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs a digital option
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsDigitalOption() 
{
    static const std::string method = "Testing Autocall vs Digital Option";
    BOOST_TEST_MESSAGE( method );

    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    fixedCoupons.back() = -1.0;

    const Schedule& callDates = fixedSchedule;
    std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    callLevels.back() = 1.0;

    std::vector<Real> callPayments( yearsToMaturity, 0.0 );
    callPayments.back() = 1.0;

    const Date strikeDate = issueDate;
    std::vector<Real> spotPrices(1, 6710.45);
    std::vector<Real> strikeLevels = spotPrices;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  QL_MAX_REAL;
    const Real redemptionUpsideGearing   =  0.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  0.0;
    const Real redemptionBarrier         =  1.0;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create a digital option
    boost::shared_ptr<CashOrNothingPayoff> payoff(new
        CashOrNothingPayoff( Option::Call, spotPrices[0], 1.0 ) );

    const Date exDate = maturityDate;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));
    EuropeanOption digitalOption( payoff, exercise );

    // Create market data
    const Size numUnderlyings = 1;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketData(
        today,
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    boost::shared_ptr<GeneralizedBlackScholesProcess> process =
        boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(
                                                      process_array->process(0));
    QL_REQUIRE(process, "Black-Scholes process required");

    digitalOption.setPricingEngine( 
        boost::shared_ptr<PricingEngine>( new AnalyticEuropeanEngine(process) ) 
        );

    const Real autocall_NPV = autocall.NPV();
    const Real option_NPV = digitalOption.NPV();
    const Real expectedResult = option_NPV/spotPrices[0];
    const Real tolerance = 1.0e-2;

    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs a call spread
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsCallSpread() 
{
    static const std::string method = "Testing Autocall vs Call Spread";
    BOOST_TEST_MESSAGE( method );

    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    fixedCoupons.back() = -2.0;

    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.0 );

    const Date strikeDate = issueDate;
    std::vector<Real> spotPrices(1, 6321.0);
    std::vector<Real> strikeLevels = spotPrices;

    const Real redemptionCash            =  2.0;
    const Real redemptionCap             =  2.6;
    const Real redemptionUpsideGearing   =  2.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  0.0;
    const Real redemptionBarrier         =  1.0;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create two options
    const Date exDate = maturityDate;
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));

    boost::shared_ptr<PlainVanillaPayoff> payoff100(
        new PlainVanillaPayoff( Option::Call, spotPrices[0] ) );
    EuropeanOption option100( payoff100, exercise );

    boost::shared_ptr<PlainVanillaPayoff> payoff130(
        new PlainVanillaPayoff( Option::Call, 1.3*spotPrices[0] ) );
    EuropeanOption option130( payoff130, exercise );

    // Create market data
    const Size numUnderlyings = 1;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketDataRealistic(
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    boost::shared_ptr<GeneralizedBlackScholesProcess> process =
        boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(
                                                      process_array->process(0));
    QL_REQUIRE(process, "Black-Scholes process required");

    option100.setPricingEngine( 
        boost::shared_ptr<PricingEngine>( new AnalyticEuropeanEngine(process) ) 
        );

    option130.setPricingEngine( 
        boost::shared_ptr<PricingEngine>( new AnalyticEuropeanEngine(process) ) 
        );

    const Real autocall_NPV = autocall.NPV();
    const Real option100_NPV = option100.NPV();
    const Real option130_NPV = option130.NPV();

    const Real expectedResult = 2.0 * (option100_NPV - option130_NPV)/spotPrices[0];
    const Real tolerance = 1.0e-2;

    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs a barrier option
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsBarrierOption() 
{
    static const std::string method = "Testing Autocall vs Barrier Option";
    BOOST_TEST_MESSAGE( method );

    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    fixedCoupons.back() = -1.0;

    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.0 );

    const Date strikeDate = issueDate;
    std::vector<Real> spotPrices(1, 6710.45);
    std::vector<Real> strikeLevels = spotPrices;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  1.0;
    const Real redemptionUpsideGearing   =  0.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  1.0;
    const Real redemptionBarrier         =  0.6;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create a barrier option
    const Date exDate = maturityDate;
    const Barrier::Type barrierType = Barrier::DownOut;
    const Real barrier = 0.6;
    const Real rebate = 0.0;
    const Option::Type type = Option::Put;
    const Real strike = spotPrices[0];
    const boost::shared_ptr<SimpleQuote> underlying( new SimpleQuote(spotPrices[0]) );

    const boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));
    const boost::shared_ptr<StrikedTypePayoff> payoff( new PlainVanillaPayoff( type, strike ) );

    BarrierOption barrierOption( barrierType, barrier, rebate, payoff, exercise );

    // Create market data
    const Size numUnderlyings = 1;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketData(
        today,
        numUnderlyings,
        riskFreeRate,
        process_array );

    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    boost::shared_ptr<GeneralizedBlackScholesProcess> process =
        boost::dynamic_pointer_cast<GeneralizedBlackScholesProcess>(
                                                      process_array->process(0));
    QL_REQUIRE(process, "Black-Scholes process required");

    barrierOption.setPricingEngine( 
        boost::shared_ptr<PricingEngine>( new AnalyticBarrierEngine(process) ) 
        );

    const Real autocall_NPV = autocall.NPV();
    const Real barrierOption_NPV = barrierOption.NPV();
    const Real expectedResult = -barrierOption_NPV/spotPrices[0];
    const Real tolerance = 1.0e-4;

    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}

// +-------+-------+-------+-------+-------+-------+-------+-------+
// | Test vs an everest option
// +-------+-------+-------+-------+-------+-------+-------+-------+
void AutocallTest::testVsEverest() 
{
    static const std::string method = "Testing Autocall vs Everest Option";
    BOOST_TEST_MESSAGE( method );
    
    const Date today(1, July, 2014);
    Date originalEvaluationDate = Settings::instance().evaluationDate();
    Settings::instance().evaluationDate() = today;

    // Create autocall
    const Integer yearsToMaturity = 6;
    const Date issueDate = today;
    const Date maturityDate = issueDate + Period( yearsToMaturity, Years );
    const std::vector<Real> notionals(1, 1.0);

    const Schedule fixedSchedule = MakeSchedule()
        .from( issueDate )
        .to( maturityDate )
        .withFrequency( Annual );

	const std::vector<Real> fixedCoupons( yearsToMaturity, 0.0 );
    const Schedule& callDates = fixedSchedule;
    const std::vector<Real> callLevels( yearsToMaturity, QL_MAX_REAL );
    const std::vector<Real> callPayments( yearsToMaturity, 0.0 );

    const Date strikeDate = issueDate;
    std::vector<Real> strikeLevels;

    const Real redemptionCash            =  1.0;
    const Real redemptionCap             =  1.0;
    const Real redemptionUpsideGearing   =  0.0;
    const Real redemptionFloor           = -1.0;
    const Real redemptionCap2            =  0.0;
    const Real redemptionDownsideGearing =  1.0;
    const Real redemptionBarrier         =  1.0;

    const Autocall::RedemptionInfo redemptionInfo(
        redemptionCash,
        redemptionCap, 
        redemptionUpsideGearing, 
        redemptionFloor,
        redemptionCap2,
        redemptionDownsideGearing,
        redemptionBarrier
        );

    Autocall autocall(
        notionals,
        fixedSchedule,
        fixedCoupons,
        callDates,
        callLevels,
        callPayments,
        strikeDate,
        strikeLevels,
        redemptionInfo
    );

    // Create an everest option 
    const Rate guarantee = 0.0;
    boost::shared_ptr<Exercise> exercise( new EuropeanExercise(maturityDate) );
    EverestOption everestOption( notionals[0], guarantee, exercise );

    // Create market data
    const Size numUnderlyings = 4;
    
    Handle<YieldTermStructure> riskFreeRate;
    boost::shared_ptr<StochasticProcessArray> process_array;

    createMarketDataRealistic(
        numUnderlyings,
        riskFreeRate,
        process_array );
    
    // Calculate NPVs
	const BigNatural seed = 86421;
	const Size fixedSamples = static_cast<Size>( std::pow(2.0,cSamplePower)-1 );
    const Size stepsPerYear = 10;

    autocall.setPricingEngine(
        MakeMCAutocallEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    everestOption.setPricingEngine(
        MakeMCEverestEngine<LowDiscrepancy>(process_array)
            .withAntitheticVariate(true)
            .withStepsPerYear(stepsPerYear)
	        .withSamples(fixedSamples)
	        .withSeed(seed)
        );

    const Real autocall_NPV = autocall.NPV();
    const Real option_NPV = everestOption.NPV();
    const Real expectedResult = option_NPV;
    const Real tolerance = 1.0e-2;
 
    Settings::instance().evaluationDate() = originalEvaluationDate;
    autocall_boost_check( method, autocall_NPV, expectedResult, tolerance );
}
