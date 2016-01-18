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

/*! \file autocall.hpp
 \brief Autocall option
 */

#ifndef quantlib_autocall_hpp
#define quantlib_autocall_hpp

#include <ql/instruments/multiassetoption.hpp>
#include <ql/instruments/basketoption.hpp>
#include <ql/time/date.hpp>
#include <ql/time/schedule.hpp>
#include <vector>

namespace QuantLib {
    
    //! Autocall
    /*! The payoff of an Autocall is computed in the following
     way:
     
     Given a group of N assets, the payout depends on either the worst
     performing of these, i.e.
     
     P(t) = min( S_k(t)/S_k(0) ) for k=0 to N - 1
     
     or on the Basket performance of these, i.e.
     
     P(t) = (1/N) * sum ( S_k(t)/S_k(0) ) for k=0 to N - 1
     
     On each coupon date :
     
     if P(t) > autocall level, then pay the autocall amount (times notional)
     and terminate, otherwise pay the fixed coupon (times notional).
     
     At maturity (if not autocalled), pay X1 + X2, where
     
     X1 = max( CashPayment, min( Cap, UpsideGearing*P(t) ) )
     X2 =
     if barrier breached then
     max( Floor, min( Cap2 , DownsideGearing*(P(T)-1) ) )
     else
     0
     
     where barrier breached == TRUE if on the maturity
     date P(T) > barrier.
     */
    
    
    class Autocall : public MultiAssetOption {
    public:
        class engine;
        class arguments;
        class results;        
        
        // Class to hold all the redemption information
        class RedemptionInfo {
        public:
            RedemptionInfo(
                           const Real& cash            =  1.0,         // Cash amount
                           const Real& cap             =  QL_MAX_REAL, // Cap
                           const Real& upsideGearing   =  1.0,         // Upside gearing
                           const Real& floor           = -1.0,         // Floor
                           const Real& cap2            =  0.0,         // 2nd cap
                           const Real& downsideGearing =  1.0,         // Downside gearing
                           const Real& barrier         =  QL_MAX_REAL  // Barrier level
            ) :
            cash_            ( cash            ),
            cap_             ( cap             ),
            upsideGearing_   ( upsideGearing   ),
            floor_           ( floor           ),
            cap2_            ( cap2            ),
            downsideGearing_ ( downsideGearing ),
            barrier_         ( barrier         )
            {}
            
            const Real& getCash()            const { return cash_;            }
            const Real& getCap()             const { return cap_;             }
            const Real& getUpsideGearing()   const { return upsideGearing_;   }
            const Real& getFloor()           const { return floor_;           }
            const Real& getCap2()            const { return cap2_;            }
            const Real& getDownsideGearing() const { return downsideGearing_; }
            const Real& getBarrier()         const { return barrier_;         }
            
        private:
            Real cash_;
            Real cap_;
            Real upsideGearing_;
            Real floor_;
            Real cap2_;
            Real downsideGearing_;
            Real barrier_;
        };
        
        Autocall(
                 const std::vector<Real>& notionals,                           // Notionals
                 const Schedule& fixedSchedule,                                // Fixed coupon schedule
                 const std::vector<Rate>& fixedCoupons,                        // Fixed coupons
                 const Schedule& callDates,                                    // Autocall dates
                 const std::vector<Real>& callLevels,                          // Autocall levels
                 const std::vector<Real>& callPayments,                        // Payments if autocalled
                 const Date& strikeDate,                                       // Strike date
                 const std::vector<Real>& strikeLevels = std::vector<Real>(),  // Strike levels
                 const RedemptionInfo& = RedemptionInfo(),                     // Redemption info
//                 const boost::shared_ptr<BasketPayoff>& underlyingType = boost::shared_ptr<BasketPayoff>(new MinBasketPayoff(boost::shared_ptr<Payoff>(new PlainVanillaPayoff(Option::Call, 0.0))) // Underlying Type
                 const boost::shared_ptr<BasketPayoff>& underlyingType = boost::shared_ptr<BasketPayoff>(new MinBasketPayoff(boost::shared_ptr<Payoff>())) // Underlying Type
        );

        void setupArguments(PricingEngine::arguments*) const;
        
        const std::vector<Real>&  getNotionals()         const { return notionals_;         }
        const Schedule&           getFixedSchedule()     const { return fixedSchedule_;     }
        const std::vector<Rate>&  getFixedCoupons()      const { return fixedCoupons_;      }
        const Schedule&           getCallDates()         const { return callDates_;         }
        const std::vector<Real>&  getCallLevels()        const { return callLevels_;        }
        const std::vector<Real>&  getCallPayments()      const { return callPayments_;      }
        const Date&               getStrikeDate()        const { return strikeDate_;        }
        const std::vector<Real>&  getStrikeLevels()      const { return strikeLevels_;      }
        const RedemptionInfo&     getRedemptionInfo()    const { return redemptionInfo_;    }
        const boost::shared_ptr<BasketPayoff>&       getUnderlyingType()    const { return underlyingType_;    }
        
    private:
        std::vector<Real>   notionals_;
        Schedule            fixedSchedule_;
        std::vector<Rate>   fixedCoupons_;
        Schedule            callDates_;
        std::vector<Real>   callLevels_;
        std::vector<Real>   callPayments_;
        Date                strikeDate_;
        std::vector<Real>   strikeLevels_;
        RedemptionInfo      redemptionInfo_;
        boost::shared_ptr<BasketPayoff>        underlyingType_;
    };
    
    class Autocall::arguments : public MultiAssetOption::arguments {
    public:
        arguments();
        void validate() const;
        
        std::vector<Real> notionals_;
        Schedule fixedSchedule_;
        std::vector<Rate> fixedCoupons_;
        Schedule callDates_;
        std::vector<Real> callLevels_;
        std::vector<Real> callPayments_;
        Date strikeDate_;
        std::vector<Real> strikeLevels_;
        RedemptionInfo redemptionInfo_;
        boost::shared_ptr<BasketPayoff> underlyingType_;
    };
    
    class Autocall::results : public MultiAssetOption::results {};
    
    class Autocall::engine
    : public GenericEngine<Autocall::arguments,
    Autocall::results> {};
}

#endif
