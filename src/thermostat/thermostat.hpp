#pragma once
#include "../output/observer.hpp"
#include "../variables.hpp"
#include "../constants.hpp"

class Thermostat {
	public:
		Variables *vars;
		Observer *obs;
		double T;
		double loop;
		double interval;
		double Tnow;
		void computeTnow(void);
		virtual void NH_zeta(void){};
		virtual void Tcontrol(int LOOP){};
		Thermostat(void){};
		~Thermostat(void){};
	private:
};

class NoseHoover : public Thermostat{
	public:
		double zeta;
		double dt;
		double Q_inv;
		void NH_zeta(void);
		void Tcontrol(int LOOP);
		NoseHoover(Variables *VAR, Observer *OBS, double DT, double temp, double Q){
			interval=100;
			loop=0;
			zeta=0;
			vars=VAR;
			obs=OBS;
			dt=DT;
			T=temp;
			Q_inv=1/Q;
		};
		~NoseHoover(void){};
	private:
};

class velocityScaling : public Thermostat {
	public:
		void NH_zeta(void){};
		void Tcontrol(int LOOP);
		velocityScaling(Variables *VAR, Observer *OBS, double temp){
			vars=VAR;
			obs=OBS;
			T=temp;
			interval=100;
			loop=0;
		};
		~velocityScaling(void){};
	private:
};

