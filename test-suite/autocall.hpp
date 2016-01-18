/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Catley Lakeman Ltd. 
 The moral rights of the author, David Rees, have been asserted.

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

#ifndef quantlib_test_autocall_hpp
#define quantlib_test_autocall_hpp

#include <boost/test/unit_test.hpp>

#include <ql/quantlib.hpp>
#include <ql/experimental/exoticoptions/autocall.hpp>

using namespace QuantLib;

class AutocallTest {
public:
    // Tests for autocall set up to replicate a fixed bond
    static void testBondBeforeIssueBeforeStrike();
    static void testBondBeforeIssueAfterStrike();
    static void testBondOnCouponDate();
    static void testBondOffCouponDate();
    static void testBondAtMaturity();
    static void testBondAfterMaturity();

    // Test for normal autocall
    static void testBeforeIssueBeforeStrike();
    static void testBeforeIssueAfterStrike();
    static void testOnCouponDate();
    static void testOffCouponDate();
    static void testAtMaturity();
    static void testAfterMaturity();

    static void testExample();

    static void testVsZeroStrikeCall();
    static void testVsCallOption();
    static void testVsDigitalOption();
    static void testVsCallSpread();
    static void testVsBarrierOption();
    static void testVsEverest();

    static boost::unit_test_framework::test_suite* suite();

private:
    // Enum used for setting the difference evaluation date tests
    enum EvaluationDateType { 
        BeforeIssueBeforeStrike,
        BeforeIssueAfterStrike,
        OnCouponDate,
        OffCouponDate,
        AtMaturity,
        AfterMaturity
    };

    // Test a autocall set up as a fixed bond
    static void testBondStyle(
        const std::string& method,
        const EvaluationDateType& evaluationDateType,
        const Real& tolerance
        );

    // Test a normal autocall 
    static void testNormal(
        const std::string& method,
        const EvaluationDateType& evaluationDateType,
        const Real& expectedResult,
        const Real& tolerance
        );

    // Create an autocall with same cashflows as fixed bond
    static void createBondStyleAutocall(
        const Date& issueDate,
        const Size& numUnderlyings,
        const bool strikeLevelsSet,
        boost::shared_ptr<Autocall>& bondStyleAutocall,
        boost::shared_ptr<FixedRateBond>& bond
        );

    // Create a normal autocall 
    static void createNormalAutocall(
        const Date& issueDate,
        const Size& numUnderlyings,
        const bool strikeLevelsSet,
        boost::shared_ptr<Autocall>& autocall
        );

    // Create market data as of valuation date
    static void createMarketData(
        const Date& valuationDate,
        const Size& numUnderlyings,
        Handle<YieldTermStructure>& riskFreeRate,
        boost::shared_ptr<StochasticProcessArray>& process_array
        );

    // Create some more realistic market data as of valuation date
    static void createMarketDataRealistic(
        const Size& numUnderlyings,
        Handle<YieldTermStructure>& riskFreeRate,
        boost::shared_ptr<StochasticProcessArray>& process_array
        );
};

#endif
