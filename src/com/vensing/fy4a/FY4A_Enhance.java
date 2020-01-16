package com.vensing.fy4a;

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.ImageIO;

import ucar.nc2.Attribute;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

/**
 * @date 2020年1月10日 上午9:01:10
 * @author vensi
 * @describe 风云四L1类数据解析
 *
 */
public class FY4A_Enhance{

	static String filePath = "data/FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20190731000000_20190731001459_4000M_V0001.HDF";
	static float Kelvin = 273.13f;

	public static void main(String[] args) throws Exception {

		NetcdfFile nc = NetcdfFile.open(filePath);

		//通道12,10.8 um
		Variable n12 = nc.findVariable("NOMChannel12");
		
		//DN值区间为[0-4096]
		short maxDN = 0;
		short minDN = 4096;
		short[] DNArray = (short[])n12.read().copyTo1DJavaArray();
		//通道12定标表,获取亮温(通过DN值索引定标表)
		Variable c12 = nc.findVariable("CALChannel12");
		float[] cal = (float[])c12.read().copyTo1DJavaArray();
		
		//求取本文件中TBB(亮温的值范围)
		for (short s : DNArray) {
			//排除无效值-1
			if(s!=-1 && s>0) {
				if(maxDN<s)maxDN = s;
				if(minDN>s)minDN = s;
			}
		}
		float maxTBB = cal[minDN];
		float minTBB = cal[maxDN];
		
		//65535无效填充值会被转换为-1
		short[][] DNs = (short[][])n12.read().copyToNDJavaArray();


		//获取起始行号及资料文件分辨率
		Attribute line_begin = nc.findGlobalAttribute("Begin_Line_Number");
		short linebegin = (short)line_begin.getValue(0);
		Attribute productName = nc.findGlobalAttribute("ProducetName");
		String product = productName.getStringValue();
		String resolution = product.split("_")[11];

		double[] geo_range = {10,54,70,140};
		double lat_S = geo_range[0], 
				lat_N = geo_range[1], 
				lon_W = geo_range[2],
				lon_E = geo_range[3];
		double step = 0.05;
		//因为循环条件取到geo_range边界值的情况，这里长度需+1
		int loncount = (int) Math.ceil((lon_E - lon_W)/step)+1;
		int latcount = (int) Math.ceil((lat_N - lat_S)/step)+1;

		//构造二维数组时，若如lonmesh[i][j]赋值，则需要初始化两个维度的长度容量
		//经纬度网格和行列号网格
		double[][] lonmesh = new double[latcount][loncount];
		double[][] latmesh = new double[latcount][loncount];
		int[][] linemesh = new int[latcount][loncount];
		int[][] columnmesh = new int[latcount][loncount];

		//dn值网格和定标亮温网格
		short[][] dnmesh = new short[latcount][loncount];
		double[][] calibratemesh = new double[latcount][loncount];
		BufferedImage bufimg=new BufferedImage(loncount, latcount, BufferedImage.TYPE_INT_ARGB);
		Graphics2D gs=bufimg.createGraphics();
		String filename  = "data/chl12-enhance.png";
		
		for(int i=0;lat_N>=lat_S;lat_N-=step,i++) {

			//每行结束后，列都从lon_W开始
			lon_W = geo_range[2];
			for(int j=0;lon_W<=lon_E;lon_W+=step,j++) {

				lonmesh[i][j] = lon_W;
				latmesh[i][j] = lat_N;

				//经纬度换算行列号
				double[] linecolumn = latlon2lc(lat_N, lon_W,resolution);
				int line =(int) Math.round(linecolumn[0])-linebegin;
				int colum =(int) Math.round(linecolumn[1]);
				linemesh[i][j] = line;
				columnmesh[i][j] = colum;

				//dn值网格和定标亮温网格赋值
				short dn = DNs[line][colum];
				dnmesh[i][j] = dn;
				calibratemesh[i][j] = cal[dn];
				
				//动画渲染，灰度值区间为后四分之三[64-255]，即[灰黑->白]，这个效果比较好
				byte g = (byte) (255 - 255*(cal[dn]-minTBB)/(maxTBB-minTBB)*0.75);
				//原始云图渲染，灰度值区间为[0-255]，即[黑->白]，这个效果黑白对比度高
				//byte g = 255 - (byte)(255*(cal[dn]-minTBB)/(maxTBB-minTBB));
				//argb
				int rgb = 255<<24;
			    rgb += g<<16;
				rgb += g<<8;
				rgb += g;
				
				//a色调增强配色方案
				rgb = enhanceColor(cal[dn],rgb);
				bufimg.setRGB(j, i, rgb);
			}
		}
		gs.dispose();  
		bufimg.flush();
		File out = new File(filename);
		System.out.println(filename);
		
		try {
			ImageIO.write(bufimg, "png", out);
		} catch (IOException e) {
			e.printStackTrace();
		}  



		System.out.println(lonmesh[0].length);
		System.out.println(lonmesh.length);
	}

	
	public static int enhanceColor(float K,int rgb) {
		//a色调增强配色方案
		int alpha = 255;
		int[][] enhance_229_224 = {
				{170,255,255},{155,220,220},{145,180,180},
				{130,150,150},{111,138,138}
		};
		int[][] enhance_224_220 = {
				{0,0,100},{0,0,150},{0,0,200},{0,0,255}
		};
		int[][] enhance_220_206 = {
				{0,100,0},{0,110,0},{0,125,0},{0,140,0},
				{0,150,0},{0,160,0},{0,170,0},{0,180,0},
				{0,190,0},{0,200,0},{0,210,0},{0,220,0},
				{0,230,0},{0,255,0}
		};
		int[][] enhance_206_199 = {
				{100,0,0},{140,0,0},{160,0,0},{180,0,0},
				{200,0,0},{220,0,0},{253,0,0}
		};
		int[][] enhance_199_191 = {
				{255,255,0},{230,230,0},{200,200,0},{190,190,0},
				{170,170,0},{150,150,0},{130,130,0},{117,117,62}
		};
		int[][] enhance_191_170 = {
				{255,255,255},{220,220,220},{200,200,200},{180,180,180},
				{160,160,160},{140,140,140},{120,120,120},{100,100,100},
				{90,90,90},{77,77,77},{0,0,50}
		};
		if(K<230&&K>225) {
			int index = (int)(229-Math.floor(K));
			if(index>4)index = 4;
			int[] color = enhance_229_224[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		
		if(K<225&&K>221) {
			int index = (int)(224-Math.floor(K));
			if(index>3)index = 3;
			int[] color = enhance_224_220[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		
		if(K<221&&K>207) {
			int index = (int)(220-Math.floor(K));
			if(index>13)index = 13;
			int[] color = enhance_220_206[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		
		if(K<207&&K>200) {
			int index = (int)(206-Math.floor(K));
			if(index>6)index = 6;
			int[] color = enhance_206_199[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		
		if(K<200&&K>191) {
			int index = (int)(199-Math.floor(K));
			if(index>7)index = 7;
			int[] color = enhance_199_191[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		if(K<191&&K>171) {
			int index = (int)(191-Math.floor(K));
			if(index>=11)index = index/2;
			int[] color = enhance_191_170[index];
			rgb = alpha<<24;
			rgb += color[0]<<16;
			rgb += color[1]<<8;
			rgb += color[2];
		}
		if(K<171) {
			rgb = alpha<<24;
			rgb += 255<<16;
			rgb += 255<<8;
			rgb += 255;
		}
		
		return rgb;
	}

	/******************************************************************
	 * 
	 * (lat, lon) → (line, column) 
	 * @warn 注意返回的是double[],非整数
	 * resolution：文件名中的分辨率{'0500M', '1000M', '2000M', '4000M'}
	 * 
	 ******************************************************************/
	public static double[] latlon2lc(double lat,double lon,String resolution) {

		double ea = 6378.137;  //地球的半长轴[km]
		double eb = 6356.7523;  //地球的短半轴[km]
		int h = 42164;  //地心到卫星质心的距离[km]
		double λD = Math.toRadians(104.7);  // 卫星星下点所在经度

		Map<String, Double> COFF = new HashMap<String, Double>();
		COFF.put("0500M",10991.5);
		COFF.put("1000M",5495.5);
		COFF.put("2000M",2747.5);
		COFF.put("4000M",1373.5);
		//列比例因子
		Map<String, Integer> CFAC = new HashMap<String, Integer>();
		CFAC.put("0500M",81865099);
		CFAC.put("1000M",40932549);
		CFAC.put("2000M",20466274);
		CFAC.put("4000M",10233137);	

		// Step1.检查地理经纬度
		// Step2.将地理经纬度的角度表示转化为弧度表示
		double lata = Math.toRadians(lat); //Math.PI*lat/180; 
		double lona = Math.toRadians(lon); //Math.PI*lon/180;

		// Step3.将地理经纬度转化成地心经纬度 
		double eb2_ea2 = Math.pow(eb, 2)/Math.pow(ea, 2);
		double λe = lona;
		double φe = Math.atan(eb2_ea2 * Math.tan(lata));

		// Step4.求Re
		double re = eb / Math.sqrt((1 - (1 - eb2_ea2) * Math.pow(Math.cos(φe),2)));

		// Step5.求r1,r2,r3
		double λe_λD = λe - λD;
		double r1 = h - re * Math.cos(φe) * Math.cos(λe_λD);
		double r2 = -re * Math.cos(φe) * Math.sin(λe_λD);
		double r3 = re * Math.sin(φe);

		// Step6.求rn,x,y
		double rn = Math.sqrt(r1*r1 + r2*r2 + r3*r3);
		double x = Math.toDegrees(Math.atan(-r2 / r1));
		double y = Math.toDegrees(Math.asin(-r3 / rn));
		// Step7.求c,l
		double c = COFF.get(resolution) + x*Math.pow(2,-16)*CFAC.get(resolution);
		double l = COFF.get(resolution) + y*Math.pow(2,-16)*CFAC.get(resolution);
		double[] linecolumn = {l,c};

		return linecolumn;
	}

	/******************************************************************
	 * 
	 * (line, column) → (lat, lon)
	 * resolution：文件名中的分辨率{'0500M', '1000M', '2000M', '4000M'}
	 * 
	 ******************************************************************/
	public static double[] lc2latlon(int line,int column,String resolution) {

		double ea = 6378.137;  //地球的半长轴[km]
		double eb = 6356.7523;  //地球的短半轴[km]
		int h = 42164;  //地心到卫星质心的距离[km]
		double λD = Math.toRadians(104.7);  // 卫星星下点所在经度

		Map<String, Double> COFF = new HashMap<String, Double>();
		COFF.put("0500M",10991.5);
		COFF.put("1000M",5495.5);
		COFF.put("2000M",2747.5);
		COFF.put("4000M",1373.5);
		//列比例因子
		Map<String, Integer> CFAC = new HashMap<String, Integer>();
		CFAC.put("0500M",81865099);
		CFAC.put("1000M",40932549);
		CFAC.put("2000M",20466274);
		CFAC.put("4000M",10233137);	

		// Step1.求x,y
		double xreg = (column-COFF.get(resolution)) / (Math.pow(2,-16)*CFAC.get(resolution));
		double yreg = (line-COFF.get(resolution)) / (Math.pow(2,-16)*CFAC.get(resolution));
		double x = Math.toRadians(xreg);
		double y = Math.toRadians(yreg);

		// Step2.求sd,sn,s1,s2,s3,sxy
		double cosx = Math.cos(x);
		double cosy = Math.cos(y);
		double siny = Math.sin(y);
		double cos2y = Math.pow(cosy, 2);
		double hcosxcosy = h * cosx * cosy;
		double cos2y_ea_eb_sin2y_2 = cos2y + Math.pow((ea/eb*siny), 2);
		double sd = Math.sqrt(Math.pow(hcosxcosy, 2) - cos2y_ea_eb_sin2y_2 * (Math.pow(h, 2)-Math.pow(ea, 2)));
		double sn = (hcosxcosy - sd)/cos2y_ea_eb_sin2y_2;
		double s1 = h - sn * cosx * cosy;
		double s2 = sn * Math.sin(x) * cosy;
		double s3 = -sn * siny;
		double sxy = Math.sqrt(s1*s1+s2*s2);

		// Step3.求lon,lat
		double lon = Math.toDegrees(Math.atan(s2/s1)+λD);
		double lat = Math.toDegrees(Math.atan(Math.pow(ea, 2)/Math.pow(eb, 2)*s3/sxy)); 
		double[] lonlat = {lon,lat};

		return lonlat;
	}

}