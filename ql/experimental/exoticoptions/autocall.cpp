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

#include <ql/experimental/exoticoptions/autocall.hpp>
#include <ql/instruments/payoffs.hpp>
#include <ql/exercise.hpp>

namespace QuantLib {
    
    Autocall::Autocall(
                       const std::vector<Real>& notionals,
                       const Schedule& fixedSchedule,
                       const std::vector<Rate>& fixedCoupons,
                       const Schedule& callDates,
                       const std::vector<Real>& callLevels,
                       const std::vector<Real>& callPayments,
                       const Date& strikeDate,
                       const std::vector<Real>& strikeLevels,
                       const Autocall::RedemptionInfo& redemptionInfo,
                       const boost::shared_ptr<BasketPayoff>& underlyingType
                       )
    : MultiAssetOption(boost::shared_ptr<Payoff>(
                                                 new PlainVanillaPayoff(Option::Call, callLevels.back())),
                       boost::shared_ptr<Exercise>(
                                                   new EuropeanExercise(callDates.endDate()))),
    notionals_      ( notionals      ),
    fixedSchedule_  ( fixedSchedule  ),
    fixedCoupons_   ( fixedCoupons   ),
    callDates_      ( callDates      ),
    callLevels_     ( callLevels     ),
    callPayments_   ( callPayments   ),
    strikeDate_     ( strikeDate     ),
    strikeLevels_   ( strikeLevels   ),
    redemptionInfo_ ( redemptionInfo ),
    underlyingType_ ( underlyingType )
    {}
    
    void Autocall::setupArguments(PricingEngine::arguments* args) const {
        MultiAssetOption::setupArguments(args);
        
        Autocall::arguments* arguments =
        dynamic_cast<Autocall::arguments*>(args);
        QL_REQUIRE(arguments != 0, "wrong argument type");
        
        arguments->notionals_       = notionals_;
        arguments->fixedSchedule_   = fixedSchedule_;
        arguments->fixedCoupons_    = fixedCoupons_;
        arguments->callDates_       = callDates_;
        arguments->callLevels_      = callLevels_;
        arguments->callPayments_    = callPayments_;
        arguments->strikeDate_      = strikeDate_;
        arguments->strikeLevels_    = strikeLevels_;
        arguments->redemptionInfo_  = redemptionInfo_;
        arguments->underlyingType_  = underlyingType_;
    }
    
    Autocall::arguments::arguments() {}
    
    void Autocall::arguments::validate() const {
        MultiAssetOption::arguments::validate();
        
        // Array sizes must be consistent
        const Size numFixed = fixedSchedule_.size();
        for( Size i=1 ; i<numFixed ; ++i )
            QL_REQUIRE( fixedSchedule_[i] > fixedSchedule_[i-1],
                       "Fixed schedule must be in ascending date order" );
        
        const Size numNotionals = notionals_.size();
        QL_REQUIRE( numNotionals==1 || numNotionals == numFixed-1,
                   "Must have either one notional or the same number of notionals and coupon dates. " <<
                   "Found " << numNotionals << " notionals and " <<
                   numFixed-1 << " coupons."
                   );
        
        const Size numFixedCoupons = fixedCoupons_.size();
        QL_REQUIRE( numFixedCoupons == 1 || numFixedCoupons == numFixed-1,
                   "Must have either one coupon or the same number of coupons and fixed coupon dates. " <<
                   "Found " << fixedCoupons_.size() << " fixed coupons and " <<
                   numFixed-1 << " fixed coupon dates."
                   );
        
        const Size numCalls = callDates_.size();
        for( Size i=1 ; i<numCalls-1 ; ++i )
            QL_REQUIRE( callDates_[i] > callDates_[i-1],
                       "Call dates must be in ascending date order" );
        
        const Size numCallLevels = callLevels_.size();
        QL_REQUIRE( numCallLevels == numCalls-1,
                   "Must have the same number of call levels and call dates. "  <<
                   "Found " << numCallLevels << " call levels and " << numCalls-1 << "call dates."
                   );
        
        const Size numCallPayments = callPayments_.size();
        QL_REQUIRE( numCallPayments == numCalls-1,
                   "Must have the same number of call payments and call dates. " <<
                   "Found " << numCallPayments << " call payments and " << numCalls-1 << "call dates."
                   );
        
        // Require that the call dates are the same as the fixed dates
        QL_REQUIRE( numFixed == numCalls,
                   "Fixed schedule and call schedule must be identical: they are different lengths" );
        
        for( Size i=0 ; i<numCalls ; ++i )
            QL_REQUIRE( callDates_[i] == fixedSchedule_[i],
                       "Fixed schedule and call schedule must be identical. " <<
                       "Found element " << i << " to be different: " <<
                       callDates_[i] << " vs " << fixedSchedule_[i]
                       );
        
        // Cannot have zero or negative strike evels
        for( Size i=0 ; i<strikeLevels_.size() ; ++i ) 
            QL_REQUIRE( strikeLevels_[i]>1.0e-6,
                       "Strike levels cannot be zero or negative. "
                       << "Invalid level found : " << strikeLevels_[i] 
                       );
    }
    
}

