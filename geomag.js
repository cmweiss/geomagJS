/*
	by Christopher Weiss (cmweiss@gmail.com)
	Suggestions for improvements are appreciated.
	
	Adapted from the geomagc software and World Magnetic Model of the NOAA
	Satellite and Information Service, National Geophysical Data Center
	http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml

	geoMagFactory() requires a world magnetic model coefficient file (available
	at the above NGDC site) which is fetched via an XMLHttpRequest and,
	therefore, requires an "Ajax" capable browser (anything this century). Due
	to web browser "same-origin" security policy, the coefficient file must be
	retrieved from the same site as the web page requesting it.

	Usage:
	geoMagFactory(url) returns a function which can compute the Earth's magnetic
	field.  If it cannot successfully fetch the URL then it returns false.
	The returned function requires two arguments, latitude and longitude (in
	decimal degrees), and, optionally, altitude in feet (default is 0), and
	time in decimal years (e.g. 2009.75 -- default is the current system time).
	
	var geoMag = geoMagFactory("http://host/dir/WMM.COF");
	var latitude = 40.0;	// decimal degrees (north is positive)
	var longitude = -80.0;	// decimal degrees (east is positive)
	var altitude = 0;	// feet (optional, default is 0)
	var time = 2006.0;	// decimal years (optional, default is the current system time)
	var myGeoMag = geoMag(latitude,longitude,altitude,time);
	var magneticVariation = myGeoMag.dec;	// Geomagnetic declination (variation) in decimal degrees -- east is positive
	var magneticDip = myGeoMag.dip;	// Geomagnetic dip in decimal degreed (down is positive)
	var magneticFieldIntensity = myGeoMag.ti;	// Total Intensity of the geomagnetic field in nanoteslas
	var magneticBH = myGeoMag.bh;	// Horizontal Intensity of the geomagnetic field in nanoteslas
	var magneticBX = myGeoMag.bx;	// North Component of the geomagnetic field in nanoteslas
	var magneticBY = myGeoMag.by;	// East Component of the geomagnetic field in nanoteslas
	var magneticBZ = myGeoMag.bz;	// Vertical Component of the geomagnetic field (down is positive)
	var lat = myGeoMag.lat;	// input latitude
	var lon = myGeoMag.lon;	//input longitude

	Differences from NGDC test values:
	Test values available at 
	http://www.ngdc.noaa.gov/geomag/WMM/data/2006TestValues_WMM2005.pdf differ
	slightly from the output of this function.  The difference seems to be due
	to rounding errors in the reporting of the test values.  NGCD's
	geomag_windows.exe reports declination in degrees and minutes rounded to 
	the nearest minute.  The output of function in this file, geomag2.js,
	converted to degrees and minutes and rounded to the nearest minute agree
	with geomag_windows.exe in all cases where the test value and output value
	differ.  It appears that the test values convert the degree & minute output
	to decimal degrees and round to the nearest hundredth of a degree.
*/

function geoMagFactory(wmm) {
	'use strict';
	function rad2deg(rad) {
		return rad * (180 / Math.PI);
	}
	function deg2rad(deg) {
		return deg * (Math.PI / 180);
	}

	var i,
		epoch = wmm.epoch,
		model = wmm.model,
		z = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
		initarray = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		maxord = 12,
		maxdeg = maxord,
		tc = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		sp = z.slice(),
		cp = z.slice(),
		pp = z.slice(),
		p = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		dp = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		a = 6378.137,
		b = 6356.7523142,
		re = 6371.2,
		a2 = a * a,
		b2 = b * b,
		c2 = a2 - b2,
		a4 = a2 * a2,
		b4 = b2 * b2,
		c4 = a4 - b4,
		c = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		cd = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		snorm = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		k = [z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice(), z.slice()],
		fn = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
		fm = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
		n, m, j, D2, flnmj;

	tc[0][0] = 0;
	sp[0] = 0.0;
	cp[0] = 1.0;
	pp[0] = 1.0;
	p[0][0] = 1;

	for (i in model) {
		if (model.hasOwnProperty(i)) {
			if (model[i].m <= model[i].n) {
				c[model[i].m][model[i].n] = model[i].gnm;
				cd[model[i].m][model[i].n] = model[i].dgnm;
				if (model[i].m !== 0) {
					c[model[i].n][model[i].m - 1] = model[i].hnm;
					cd[model[i].n][model[i].m - 1] = model[i].dhnm;
				}
			}
		}
	}
	wmm = null;
	model = null;

	/* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
	snorm[0][0] = 1;

	for (n = 1; n <= maxord; n++) {
		snorm[0][n] = snorm[0][n - 1] * (2 * n - 1) / n;
		j = 2;
		m = 0;
		D2 = (n - m + 1);
		while (D2 > 0) { 
//			for (m = 0, var D2 = (n - m + 1); D2 > 0; D2--, m++) {
			k[m][n] = (((n - 1) * (n - 1)) - (m * m)) / ((2 * n - 1) * (2 * n - 3));
			if (m > 0) {
				flnmj = ((n - m + 1) * j) / (n + m);
				snorm[m][n] = snorm[m - 1][n] * Math.sqrt(flnmj);
				j = 1;
				c[n][m - 1] = snorm[m][n] * c[n][m - 1];
				cd[n][m - 1] = snorm[m][n] * cd[n][m - 1];
			}
			c[m][n] = snorm[m][n] * c[m][n];
			cd[m][n] = snorm[m][n] * cd[m][n];
			D2--; m++;
		}
//			}
	}
	k[1][1] = 0.0;

	return function (dlat, dlon, h, time) {
		var now,
			alt = h / 3280.8399, // convert h (in feet) to meters
			otime = -1000.0,
			oalt = otime,
			olat = otime,
			olon = otime,
			dt = time - epoch,
			glat = dlat,
			glon = dlon,
			rlat = deg2rad(glat),
			rlon = deg2rad(glon),
			srlon = Math.sin(rlon),
			srlat = Math.sin(rlat),
			crlon = Math.cos(rlon),
			crlat = Math.cos(rlat),
			srlat2 = srlat * srlat,
			crlat2 = crlat * crlat,
			q, q1, q2, ct, st,
			r2, r, d, ca, sa, aor, ar, br, bt, bp, bpp, par, temp1, temp2, parp, D4,
			bx, by, bz, bh, ti, dec, dip, gv;
		if (h === undefined) { h = 0; }
		if (time === undefined) {
			now = new Date();
			time = now.getFullYear() + (now.getMonth() + ((now.getDate() - 1) / 31)) / 12;
		}

		sp[1] = srlon;
		cp[1] = crlon;

		/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
		if (alt !== oalt || glat !== olat) {
			q = Math.sqrt(a2 - c2 * srlat2);
			q1 = alt * q;
			q2 = ((q1 + a2) / (q1 + b2)) * ((q1 + a2) / (q1 + b2));
			ct = srlat / Math.sqrt(q2 * crlat2 + srlat2);
			st = Math.sqrt(1.0 - (ct * ct));
			r2 = (alt * alt) + 2.0 * q1 + (a4 - c4 * srlat2) / (q * q);
			r = Math.sqrt(r2);
			d = Math.sqrt(a2 * crlat2 + b2 * srlat2);
			ca = (alt + d) / r;
			sa = c2 * crlat * srlat / (r * d);
		}
		if (glon !== olon) {
			for (m = 2; m <= maxord; m++) {
				sp[m] = sp[1] * cp[m - 1] + cp[1] * sp[m - 1];
				cp[m] = cp[1] * cp[m - 1] - sp[1] * sp[m - 1];
			}
		}
		aor = re / r;
		ar = aor * aor;
		br = 0.0;
		bt = 0.0; 
		bp = 0.0;
		bpp = 0.0;

		for (n = 1; n <= maxord; n++) {
			ar = ar * aor;
			m = 0; D4 = (n + m + 1); while (D4 > 0) {

		/*
				COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
				AND DERIVATIVES VIA RECURSION RELATIONS
		*/
				if (alt !== oalt || glat !== olat) {
					if (n === m) {
						p[m][n] = st * p[m - 1][n - 1];
						dp[m][n] = st * dp[m - 1][n - 1] + ct * p[m - 1][n - 1];
					} else if (n === 1 && m === 0) {
						p[m][n] = ct * p[m][n - 1];
						dp[m][n] = ct * dp[m][n - 1] - st * p[m][n - 1];
					} else if (n > 1 && n !== m) {
						if (m > n - 2) { p[m][n - 2] = 0; }
						if (m > n - 2) { dp[m][n - 2] = 0.0; }
						p[m][n] = ct * p[m][n - 1] - k[m][n] * p[m][n - 2];
						dp[m][n] = ct * dp[m][n - 1] - st * p[m][n - 1] - k[m][n] * dp[m][n - 2];
					}
				}
		/*
				TIME ADJUST THE GAUSS COEFFICIENTS
		*/

				tc[m][n] = c[m][n] + dt * cd[m][n];
				if (m !== 0) { tc[n][m - 1] = c[n][m - 1] + dt * cd[n][m - 1]; }

		/*
				ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
		*/
				par = ar * p[m][n];
				if (m === 0) {
					temp1 = tc[m][n] * cp[m];
					temp2 = tc[m][n] * sp[m];
				} else {
					temp1 = tc[m][n] * cp[m] + tc[n][m - 1] * sp[m];
					temp2 = tc[m][n] * sp[m] - tc[n][m - 1] * cp[m];
				}
				bt = bt - ar * temp1 * dp[m][n];
				bp += (fm[m] * temp2 * par);
				br += (fn[n] * temp1 * par);
		/*
					SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
		*/
				if (st === 0.0 && m === 1) {
					if (n === 1) { pp[n] = pp[n - 1]; } else { pp[n] = ct * pp[n - 1] - k[m][n] * pp[n - 2]; }
					parp = ar * pp[n];
					bpp += (fm[m] * temp2 * parp);
				}
				D4--; m++;
			}
		}

		if (st === 0.0) { bp = bpp; } else { bp /= st; }
		/*
			ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
			GEODETIC COORDINATES
		*/
		bx = -bt * ca - br * sa;
		by = bp;
		bz = bt * sa - br * ca;

		/*
			COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
			TOTAL INTENSITY (TI)
		*/
		bh = Math.sqrt((bx * bx) + (by * by));
		ti = Math.sqrt((bh * bh) + (bz * bz));
		dec = rad2deg(Math.atan2(by, bx));
		dip = rad2deg(Math.atan2(bz, bh));
		/*
			COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
			GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
			(I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)
			OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
		*/
		gv = -999.0;
		if (Math.abs(glat) >= 55.0) {
			if (glat > 0.0 && glon >= 0.0) { gv = dec - glon; }
			if (glat > 0.0 && glon < 0.0) { gv = dec + Math.abs(glon); }
			if (glat < 0.0 && glon >= 0.0) { gv = dec + glon; }
			if (glat < 0.0 && glon < 0.0) { gv = dec - Math.abs(glon); }
			if (gv > +180.0) { gv -= 360.0; }
			if (gv < -180.0) { gv += 360.0; }
		}
//		otime = time;
//		oalt = alt;
//		olat = glat;
//		olon = glon;

		return {dec: dec, dip: dip, ti: ti, bh: bh, bx: bx, by: by, bz: bz, lat: dlat, lon: dlon};
	};
}