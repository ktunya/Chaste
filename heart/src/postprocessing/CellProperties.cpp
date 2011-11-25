/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include <cmath>
#include <sstream>
#include <cassert>

#include "CellProperties.hpp"
#include "Exception.hpp"


enum APPhases { BELOWTHRESHOLD , ABOVETHRESHOLD };

void CellProperties::CalculateProperties()
{
    // Check we have some suitable data to process
    if (mrTime.size() < 1)
    {
        EXCEPTION("Insufficient time steps to calculate physiological properties.");
    }

    if (mrTime.size() != mrVoltage.size())
    {
        EXCEPTION("Time and Voltage series should be the same length. Time.size() = " << mrTime.size() << ", Voltage.size() = " << mrVoltage.size());
    }


    double max_upstroke_velocity = -DBL_MAX;
    double current_time_of_upstroke_velocity = 0;
    double current_resting_value=DBL_MAX;
    double current_peak=-DBL_MAX;
    double current_minimum_velocity=DBL_MAX;
    double prev_voltage_derivative=0;
    unsigned ap_counter = 0;
    unsigned counter_of_plateau_depolarisations = 0;
    //boolean to keep track whether we are switching phase from BELOWTHRESHOLD to ABOVETHRESHOLD
    bool switching_phase = false;
    bool found_a_flat_bit=false;
    APPhases ap_phase = BELOWTHRESHOLD;

    unsigned time_steps = mrTime.size()-1; //The number of time steps is the number of intervals

    double v = mrVoltage[0];
    double t = mrTime[0];
    double prev_v = v;
    double prev_t = t;
    double voltage_derivative;
    const double resting_potential_gradient_threshold = 1e-2; /// \todo #1495 a horrible magic number but seems to work OK.

    for (unsigned i=1; i<=time_steps; i++)
    {
        v = mrVoltage[i];
        t = mrTime[i];
        voltage_derivative = (t == prev_t) ? 0.0 : (v - prev_v) / (t - prev_t);

        // Look for the max upstroke velocity and when it happens (could be below or above threshold).
        if (voltage_derivative >= max_upstroke_velocity)
        {
            max_upstroke_velocity = voltage_derivative;
            current_time_of_upstroke_velocity = t;

        }

        switch (ap_phase)
        {
            case BELOWTHRESHOLD:
                //while below threshold, find the resting value by checking where the velocity is minimal
                //i.e. when it is flattest. If we can't find a flat bit, instead go for the minimum voltage
                //seen before the threshold.
                if (fabs(voltage_derivative)<=current_minimum_velocity && fabs(voltage_derivative)<=resting_potential_gradient_threshold)
                {
                    current_minimum_velocity=fabs(voltage_derivative);
                    current_resting_value = prev_v;
                    found_a_flat_bit=true;
                }
                else if(prev_v < current_resting_value && !found_a_flat_bit)
                {
                    current_resting_value = prev_v;
                }

                // If we cross the threshold, this counts as an AP
                if ( v>mThreshold && prev_v <= mThreshold )
                {
                    //register the resting value and re-initialise the minimum velocity
                    mRestingValues.push_back(current_resting_value);
                    current_minimum_velocity = DBL_MAX;
                    current_resting_value = DBL_MAX;

                    //Register the onset time. Linear interpolation.
                    mOnsets.push_back(prev_t + (t-prev_t)/(v-prev_v)*(mThreshold-prev_v));

                    //If it is not the first AP, calculate cycle length for the last two APs
                    if (ap_counter>0)
                    {
                        mCycleLengths.push_back( mOnsets[ap_counter]-mOnsets[ap_counter-1] );
                    }

                    switching_phase = true;
                    found_a_flat_bit = false;
                    ap_phase = ABOVETHRESHOLD;
                    // no break here - deliberate fall through to next case
                }
                else
                {
                    break;
                }

            case ABOVETHRESHOLD:
                //While above threshold, look for the peak potential for the current AP
                if (v>current_peak)
                {
                   current_peak = v;
                }

                // we check whether we have above threshold depolarisations
                // and only if if we haven't just switched from below threshold at this time step.
                // The latter is to avoid recording things depending on resting behaviour (in case of sudden upstroke from rest)
                if (prev_voltage_derivative<=0 && voltage_derivative>0 && !switching_phase)
                {
                    counter_of_plateau_depolarisations++;
                }

                // From the next time step, we are not "switching phase" any longer
                // (we want to check for above threshold deolarisations)
                switching_phase = false;

                // If we cross the threshold again, the AP is over
                // and we register all the parameters.
                if ( v<mThreshold && prev_v >= mThreshold )
                {
                    //register peak value for this AP
                    mPeakValues.push_back(current_peak);
                    //Re-initialise the current_peak.
                    current_peak = mThreshold;

                    //register maximum upstroke velocity for this AP
                    mMaxUpstrokeVelocities.push_back(max_upstroke_velocity);
                    //re-initialise max_upstroke_velocity
                    max_upstroke_velocity = -DBL_MAX;

                    //register time when maximum upstroke velocity occurred for this AP
                    mTimesAtMaxUpstrokeVelocity.push_back(current_time_of_upstroke_velocity);
                    //re-initialise current_time_of_upstroke_velocity=t;
                    current_time_of_upstroke_velocity = 0.0;

                    mCounterOfPlateauDepolarisations.push_back(counter_of_plateau_depolarisations);

                    //update the counters.
                    ap_counter++;
                    ap_phase = BELOWTHRESHOLD;

                    //reinitialise counter of plateau depolarisations
                    counter_of_plateau_depolarisations = 0;
                }
                break;
        }
        prev_v = v;
        prev_t = t;
        prev_voltage_derivative = voltage_derivative;
    }


    // One last check. If the simulation ends halfway through an AP
    // i.e. if the vectors of onsets has more elements than the vectors
    // of peak and upstroke properties (that are updated at the end of the AP),
    // then we register the peak and upstroke values so far
    // for the last incomplete AP.
    if (mOnsets.size()>mMaxUpstrokeVelocities.size())
    {
        mMaxUpstrokeVelocities.push_back(max_upstroke_velocity);
        mPeakValues.push_back(current_peak);
        mTimesAtMaxUpstrokeVelocity.push_back(current_time_of_upstroke_velocity);
        mUnfinishedActionPotentials = true;
    }
}


std::vector<double> CellProperties::CalculateActionPotentialDurations(const double percentage)
{
    CheckExceededThreshold();

    double prev_v = mrVoltage[0];
    unsigned APcounter=0;//will keep count of the APDs that we calculate
    bool apd_is_calculated=true;//this will ensure we hit the target only once per AP.
    std::vector<double> apds;
    double target = DBL_MAX;
    bool apd_starting_time_found=false;
    double apd_start_time=DBL_MAX;

    double t;
    double v;
    for (unsigned i=1; i<mrTime.size(); i++)
    {
        t = mrTime[i];
        v = mrVoltage[i];

        //First we make sure we stop calculating after the last AP has been calculated
        if (APcounter<mPeakValues.size())
        {
            //Set the target potential
            target = mRestingValues[APcounter]+0.01*(100-percentage)*(mPeakValues[APcounter]-mRestingValues[APcounter]);

            //if we reach the peak, we need to start to calculate an APD
            if (fabs(v-mPeakValues[APcounter])<=1e-6)
            {
                apd_is_calculated = false;
            }

            // Start the timing where we first cross the target voltage
            if ( prev_v<v && prev_v<=target && v>=target && apd_starting_time_found==false)
            {
                // Linear interpolation of target crossing time.
                apd_start_time=t+( (target-prev_v)/(v-prev_v) )*(t-mrTime[i-1]);
                apd_starting_time_found = true;
            }

            //if we hit the target while repolarising
            //and we are told this apd is not calculated yet.
            if ( prev_v>v && prev_v>=target && v<=target && apd_is_calculated==false)
            {
                // Linear interpolation of target crossing time.
                apds.push_back (t - apd_start_time + ( (target-prev_v)/(v-prev_v) )*(t-mrTime[i-1]) );
                APcounter++;
                apd_is_calculated = true;
                apd_starting_time_found = false;
            }
        }
        prev_v = v;
    }
    if (apds.size() == 0)
    {
        EXCEPTION("No full action potential was recorded");
    }
    return apds;
}

std::vector<double> CellProperties::GetActionPotentialAmplitudes()
{
    CheckExceededThreshold();
    unsigned size = mPeakValues.size();
    std::vector<double> amplitudes(size);
    for (unsigned i=0; i< size ;i++)
    {
        amplitudes[i] = (mPeakValues[i] - mRestingValues[i]);
    }
    return amplitudes;
}

void CellProperties::CheckExceededThreshold(void)
{
    // mOnsets and mRestingValues are all
    // set at the same time, so checking one should suffice.
    if (mOnsets.empty())
    {
        // possible false error here if the simulation started at time < 0
        EXCEPTION("AP did not occur, never exceeded threshold voltage.");
    }
}

void CellProperties::CheckReturnedToThreshold(void)
{
    // mPeakValues, mMaxUpstrokeVelocities and mTimesAtMaxUpstrokeVelocity are all
    // set at the same time so checking one should suffice.
    if (mMaxUpstrokeVelocities.empty())
    {
        EXCEPTION("AP did not occur, never descended past threshold voltage.");
    }
}

//
// The Get std::vector<double> methods
//

std::vector<double> CellProperties::GetCycleLengths()
{
    return mCycleLengths;
}

std::vector<double> CellProperties::GetRestingPotentials()
{
    CheckExceededThreshold();
    return mRestingValues;
}

std::vector<double> CellProperties::GetPeakPotentials()
{
    CheckReturnedToThreshold();
    return mPeakValues;
}

std::vector<double> CellProperties::GetMaxUpstrokeVelocities()
{
    CheckReturnedToThreshold();
    return mMaxUpstrokeVelocities;
}

std::vector<double> CellProperties::GetTimesAtMaxUpstrokeVelocity()
{
    CheckReturnedToThreshold();
    return mTimesAtMaxUpstrokeVelocity;
}

std::vector<double> CellProperties::GetAllActionPotentialDurations(const double percentage)
{
    return CalculateActionPotentialDurations(percentage);
}

std::vector<unsigned> CellProperties::GetNumberOfAboveThresholdDepolarisationsForAllAps()
{
    return mCounterOfPlateauDepolarisations;
}
//
// The Get <double> methods
//

double CellProperties::GetLastPeakPotential()
{
    CheckReturnedToThreshold();
    return mPeakValues.back();
}

double CellProperties::GetLastCompletePeakPotential()
{
    CheckReturnedToThreshold();
    double peak_value;
    if (mUnfinishedActionPotentials && mPeakValues.size()>1u)
    {
        peak_value = mPeakValues[mPeakValues.size()-2u];
    }
    else if  (mUnfinishedActionPotentials && mPeakValues.size()==1u)
    {
        EXCEPTION("No peak potential matching a full action potential was recorded.");
    }
    else
    {
        peak_value = mPeakValues.back();
    }
    return peak_value;
}

double CellProperties::GetLastMaxUpstrokeVelocity()
{
    CheckReturnedToThreshold();
    return mMaxUpstrokeVelocities.back();
}

double CellProperties::GetLastCompleteMaxUpstrokeVelocity()
{
    CheckReturnedToThreshold();
    double max_upstroke;
    if (mUnfinishedActionPotentials && mMaxUpstrokeVelocities.size()>1u)
    {
        max_upstroke = mMaxUpstrokeVelocities[mMaxUpstrokeVelocities.size()-2u];
    }
    else if  (mUnfinishedActionPotentials && mMaxUpstrokeVelocities.size()==1u)
    {
        EXCEPTION("No MaxUpstrokeVelocity matching a full action potential was recorded.");
    }
    else
    {
        max_upstroke = mMaxUpstrokeVelocities.back();
    }
    return max_upstroke;
}

double CellProperties::GetTimeAtLastMaxUpstrokeVelocity()
{
    CheckReturnedToThreshold();
    return mTimesAtMaxUpstrokeVelocity.back();
}

double CellProperties::GetTimeAtLastCompleteMaxUpstrokeVelocity()
{
    CheckReturnedToThreshold();
    double max_upstroke_time;
    if (mUnfinishedActionPotentials && mTimesAtMaxUpstrokeVelocity.size()>1u)
    {
        max_upstroke_time = mTimesAtMaxUpstrokeVelocity[mTimesAtMaxUpstrokeVelocity.size()-2u];
    }
    else if  (mUnfinishedActionPotentials && mTimesAtMaxUpstrokeVelocity.size()==1u)
    {
        EXCEPTION("No TimeAtMaxUpstrokeVelocity matching a full action potential was recorded.");
    }
    else
    {
        max_upstroke_time = mTimesAtMaxUpstrokeVelocity.back();
    }
    return max_upstroke_time;
}

double CellProperties::GetLastActionPotentialDuration(const double percentage)
{
    std::vector<double> apds = CalculateActionPotentialDurations(percentage);
    // We tested this returns a non-empty vector in the method.
    return apds.back();
}

unsigned CellProperties::GetNumberOfAboveThresholdDepolarisationsForLastAp()
{
    return mCounterOfPlateauDepolarisations.back();
}

