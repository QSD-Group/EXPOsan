// Rebuilds the researched snapshot; does not crawl or refresh websites.
import fs from 'node:fs/promises';
import path from 'node:path';
import os from 'node:os';
import { fileURLToPath } from 'node:url';
import { createRequire } from 'node:module';
const modules=process.env.ARTIFACT_TOOL_NODE_MODULES || path.join(os.homedir(),'.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules');
const require=createRequire(path.join(modules,'_resolver.cjs'));
const {Workbook,SpreadsheetFile}=require('@oai/artifact-tool');
const out=path.dirname(fileURLToPath(import.meta.url));
const previews=path.resolve(out,'../../outputs/chemical-prices-20261031');
const date=new Date('2026-10-31T12:00:00Z');
const u={
 fengda:'https://fengdachem.en.made-in-china.com/product/bFAfGwOvlUpd/China-Maleic-Acid-Price-99-Min-White-Crystal-110-16-7-C4h4o4-Maleic-Acid.html',
 jinao:'https://jinao2025.en.made-in-china.com/product/ITpRAQnMDrYi/China-Manufacturer-Supply-Bulk-25kg-Maleic-Acid-CAS-110-16-7-Maleic-Acid-Price.html',
 dobo:'https://dobo2021.en.made-in-china.com/product/tGoYIXVMhqhp/China-Manufacturer-Supply-Bulk-25kg-Maleic-Acid-CAS-110-16-7-Maleic-Acid-Price.html',
 acidUnits:'https://www.made-in-china.com/products-search/hot-china-products/Acid_Maleic.html',
 longze:'https://longzechem66.en.made-in-china.com/product/kQMrsfavVzpn/China-Food-Grade-Maleic-Acid-99-CAS-No-110-16-7.html',
 xima:'https://33597786e6bf458b.en.made-in-china.com/product/TaoUMJqcvVRz/China-99-Maleic-Acid-CAS-110-16-7.html',
 eastCosmetic:'https://dl-eastland.en.made-in-china.com/product/xnzUiVJwYCrF/China-Cosmetic-Raw-Materials-Chemicals-High-Quality-Ammonium-Lactate-CAS-515-98-0.html',
 eastFactory:'https://dl-eastland.en.made-in-china.com/product/dmKYQRIDglUw/China-Factory-Supply-Ammonium-Lactate-CAS-515-98-0.html',
 kejian:'https://www.alibaba.com/product-detail/60-Ammonium-Lactate-Solution-CAS-515_1601536283370.html',
 sinsur:'https://sinsurfoodadditive.en.made-in-china.com/product/SByxQhkOlqRv/China-China-Manufacturer-High-Purity-Food-Additive-Ammonium-Citrate.html',
 kolod:'https://jskolod.en.made-in-china.com/product/xsXnCGQUYpkZ/China-High-Purity-Ammonium-Citrate-Food-Grade-Bp.html',
 ditai:'https://ditaichem.en.made-in-china.com/product/uZeJWiXUZpko/China-High-Quality-Ammonium-Citrate-CAS-3458-72-8.html',
 essence:'https://www.alibaba.com/product-detail/Ammonium-Citrate-99-Food-Grade-CAS_1601785980499.html',
 eastCitrate:'https://dl-eastland.en.made-in-china.com/product/oElUrwhxbBpZ/China-Diammonium-Citrate-CAS-No-3012-65-5-Ammonium-Citrate-Dibasic.html',
 randa:'https://hbranda.en.made-in-china.com/product/FrPRTLQAuEkX/China-High-Quality-Ammonium-Citrate-Tribasic-for-Efficient-Industrial-Applications-CAS-3458-72-8.html'
};

// One row per supplier product, retaining only its lowest advertised USD quote.
Object.assign(u,{
 surest:'https://surest.en.made-in-china.com/product/StOUohQbERVZ/China-99-Maleic-Acid-CAS-110-16-7-Maleate-White-Powder.html',
 sinobio:'https://sinobioch.en.made-in-china.com/product/oZSthcLvEnAC/China-High-Purity-99-Maleic-Acid-CAS-110-16-7.html',
 mikemLactate:'https://www.alibaba.com/product-detail/High-Purity-Ammonium-Lactate-CAS-515_1601605987113.html',
 dayangLactate:'https://www.lookchem.com/product_2-Hydroxypropanoic-acid-monoammonium-salt/15111346.html',
 xinyang:'https://www.ec21.com/product-details/Ammonium-Hydrogen-Citrate-Diammonium-Citrate--10979734.html'
});
const rows=[];
function add(id,chemical,form,cas,assay,price,unit,moq,moqu,tier,pack,supplier,verification,flag,notes,url,support=''){
 rows.push([id,chemical,form,cas,assay,price,'USD',unit,moq,moqu,tier,pack,supplier,verification,date,'Pending',flag,notes,url,support]);
}
add('MA01','Maleic acid','Free acid crystals','110-16-7','Total acid ≥99%',2350,'Ton',20,'Tons','Lower bound; quantity break not specified','25 kg bag','Mudanjiang Fengda Chemicals Imp. & Exp. Corp.','Opened supplier product page','None identified','Price denominator verified against matching marketplace offer; lower bound retained.',u.fengda,u.acidUnits);
add('MA03','Maleic acid','Free acid powder','110-16-7','Specification 99%; quality ≥99.5%; assay 99.0–101.0%',2.2,'kg',1000,'kg','Lower bound; quantity break not specified','25 kg bag / 200 kg option','NINGBO DONGBO NEW ENERGY CO., LTD.','Opened supplier product page','Different assay fields','Assay claims differ; confirm specification before selection. Packaging options are alternatives. Unit verified in matching marketplace offer.',u.dobo,u.acidUnits);
add('MA05','Maleic acid','Free acid powder','110-16-7','C4H4O4 assay 99.0–101.0 wt.%',2100,'Ton',1,'Ton','Lower bound; quantity break not specified','25 kg bag','Shandong Longze Chemical Co.,Ltd','Opened supplier product page in browser','Analytical assay exceeds 100%','Seller analytical assay as published; not interpreted as a physical mass fraction. Unit verified in matching offer.',u.longze,u.acidUnits);
add('MA07','Maleic acid','Free acid crystals','110-16-7','≥99%',2050,'Ton',1,'Ton','100+ Tons required for this price','25 kg woven / composite bag','Shandong Xima Supply Chain Management Co., Ltd.','Opened supplier product page in browser','None identified','Lowest published tier only; overall MOQ is 1 Ton but selected price needs 100+ Tons. MA06 consolidated into MA07.',u.xima);
add('MA08','Maleic acid','Free acid white powder','110-16-7','99% minimum',5,'kg',25,'kg','Fixed quote; ocean-order MOQ 25 kg','25 kg bag','Ningbo Surest Chemical Co., Ltd.','Opened supplier product page in browser','Density unit not stated','USD amount, purity and packaging on product page; /kg denominator verified in exact-title marketplace offer.',u.surest,u.acidUnits);
add('MA09','Maleic acid','Free acid powder','110-16-7','99% (product title)',2,'kg',1000,'kg','Lower bound; quantity break not specified','25 kg bag / 500 kg bag','Dalian Sinobio Chemistry Co., Ltd.','Opened supplier product page in browser','None identified','Alternative bag sizes; /kg denominator verified in exact-title marketplace offer.',u.sinobio,u.acidUnits);
add('AL08','Ammonium lactate','Liquid; colorless / light yellow','515-98-0','99% purity; assay basis not specified',7,'kg',25,'kg','≥750 kg required for this price','25 kg ordering batch; container not specified','Mikem Chemical (Shandong) Co., Ltd.','Opened Alibaba supplier product page in browser','Purity basis not specified','Kilogram-priced quantity tiers verified visually. Lowest tier retained; batch size 25 kg. Seller calls this purity, not active content; do not assume dry-solids concentration.',u.mikemLactate);
add('AL09','Ammonium lactate','2-Hydroxypropanoic acid monoammonium salt','515-98-0','98% minimum',3,'Kilogram',1,'Kilogram','Fixed quote; bulk quantity break not stated','According to client requirements; size not specified','Dayang Chem (Hangzhou) Co.,Ltd.','Opened supplier product page in browser','Physical form / assay basis unspecified','Explicit USD/Kilogram and purity. Industrial raw material offer with stated 1–100 metric ton/month capacity; small MOQ is not a bulk-discount guarantee.',u.dayangLactate);
add('AC12','Diammonium citrate','Dibasic citrate; powder / granules','3012-65-5','99% purity',2.2,'Kilogram',100,'Kilogram','≥10,000 kg required for this price','25 kg; container type not specified','Ningxiang Xinyang Chemical','Indexed text of supplier product page; live browser blocked','Live page needs recheck','Primary EC21 offer text retrieved by search includes CAS, purity, explicit US$/Kilogram, packaging and tiers. Live Cloudflare verification prevented direct browser access. Template condition says Used–Very Good; supplier confirmation recommended.',u.xinyang);

const extraURLs={
 youze:'https://ouze.lookchem.com/products/CasNo-515-98-0-Ammonium-lactate-CAS-515-98-0-the-cheapest-price-28349813.html',
 sankai:'https://sankai.lookchem.com/products/CasNo-515-98-0-2-Hydroxypropanoic-acid-monoammonium-salt-CAS515-98-0-28715117.html',
 dingyan:'https://www.lookchem.com/product_High-purity-Diammonium-citrate-Ammonium-Citrate-Dibasic-with-high-quality-and-best-price-cas-3012-65-5/20512304.html',
 dayangCitrate:'https://www.lookchem.com/product_Ammonium-citrate-dibasic-Manufacturer-High-quality-Best-price-In-stock/2048955.html',
 gihi:'https://gihichem.lookchem.com/products/CasNo-3012-65-5-Ammonium-citrate-dibasic-32555664.html'
};
add('AL10','Ammonium lactate','White crystalline powder','515-98-0','99% purity; analytical basis unspecified',9,'Kilogram',1,'Kilogram','Lower bound; quantity break not specified','Customized; aluminum foil bag / fluorinated bottle / other options; mass not stated','Hebei Youze Biotechnology Co., Ltd.','Opened supplier product page in browser','Supplier page marked expired','Page explicitly states USD/Kilogram and product purity. Supplier display name includes expird/expired; current availability and quotation date need confirmation. No mass or volume inferred for customized packaging.',extraURLs.youze);
add('AL11','Ammonium lactate','White powder','515-98-0','98% product purity; basis unspecified',289,'Kilogram',1,'Kilogram','101–1,000 kg required for this price','Fiber can; net mass not stated','Hebei Sankai Chemical Technology Co., Ltd.','Opened supplier product page in browser','High price outlier; assay wording differs','Explicit USD/Kilogram tier for 101–1,000 kg. Much higher than other lactate offers; keep separate pending selection. Product field states 98%; generic marketing says all above 98.5%. Density copied-looking / basis unresolved; retained as seller claim.',extraURLs.sankai);
add('AC13','Diammonium citrate','Colorless crystals / crystalline granules','3012-65-5','Assay ≥98%; seller test result 98.52%',2,'Kilogram',1,'Kilogram','Lower bound; quantity break not specified','25 kg net paper drum; smaller-order bag/bottle options','Hangzhou Dingyan Chem Co., Ltd.','Opened supplier product page in browser','Assay basis not stated','CAS and molecular formula consistent. Product-specific 25 kg paper drum agrees with generic ≥25 kg drum/bag options. Lower quote bound only; seller does not tie it to a specific order size.',extraURLs.dingyan);
add('AC14','Diammonium citrate','Ammonium citrate dibasic','3012-65-5','98% purity',3,'Kilogram',1,'Kilogram','Fixed quote; bulk quantity break not stated','According to client requirements; size not specified','Dayang Chem (Hangzhou) Co.,Ltd.','Opened supplier product page in browser','Physical form / assay basis unspecified','Explicit USD/Kilogram industrial raw-material offer, with stated 1–100 metric ton/month capacity. CAS and formula C6H8O7·2H3N consistent; custom packaging without conflicting claims.',extraURLs.dayangCitrate);
add('AC15','Diammonium citrate','Powder; ammonium citrate dibasic','3012-65-5','99% purity',15,'Kilogram',1,'Kilogram','1–1,000 kg quoted at this price','25 kg/drum or customized; smaller-order packaging options','GIHI CHEMICALS CO.,LIMITED','Opened supplier product page in browser','Assay basis not stated','USD/Kilogram explicitly shown for 1–1,000 kg. Product details specify 25 kg/drum or customized. Generic packaging section offers different containers by quantity; treated as alternatives.',extraURLs.gihi);

const properties=[];
function prop(ids,name,lo,hi,unit,conditions,limitation,url){for(const id of ids)properties.push([id,name,lo,hi,unit,conditions,limitation,url,date]);}
prop(['MA01'],'Molecular weight',116.07,116.07,'g/mol','Free acid','Seller-reported',u.fengda);
prop(['MA01'],'Water solubility',79,79,'g/100 g','25 °C; denominator basis unclear','Seller-reported',u.fengda);
prop(['MA03'],'Density',1.59,1.59,'g/mL','25 °C; seller cites literature','Seller-reported',u.dobo);
prop(['MA03'],'Melting point',137,140,'°C','As listed','Seller-reported; confirm before use',u.dobo);
prop(['MA07'],'Density',1.59,1.59,'g/cm³','20 °C','Seller-reported',u.xima);
prop(['MA07'],'Melting point',130,131,'°C','As listed','Seller-reported',u.xima);
prop(['MA08'],'Density',1.59,1.59,'Not stated','Temperature not stated','Do not assign or infer units',u.surest);
prop(['MA08'],'Melting point',134,138,'°C','As listed','Seller-reported',u.surest);
prop(['MA08'],'Water solubility',790,790,'g/L','25 °C','Seller-reported',u.surest);
prop(['MA08'],'Molecular weight',116.07,116.07,'g/mol','As listed','Seller-reported',u.surest);

for(const id of ['AL10','AL11']){
 const source=id==='AL10'?extraURLs.youze:extraURLs.sankai;
 prop([id],'Density',1.054,1.054,'g/mL','25 °C; seller also describes powder','Physical-form / concentration basis unresolved; do not use for conversion',source);
 prop([id],'Melting point',91,94,'°C','As listed','Seller-reported; not independently validated',source);
 prop([id],'Molecular weight',107.108,107.108,'g/mol','Ammonium lactate','Seller-reported',source);
}

const log=[];
function reject(ids,chemical,reason,url){log.push([ids,chemical,'Excluded / not priced',reason,url,date]);}
reject('MA02','Maleic acid','Title says 25 kg, while packaging attributes give drums / IBC; excluded for inconsistent packaging.',u.jinao);
reject('MA06','Maleic acid','Higher Xima quantity tier consolidated into MA07; one lowest price per supplier product.',u.xima);
reject('AL01 / AL07','Ammonium lactate','Only generic active-substance content, not product purity; price denominator unclear; AL07 also bag/solution conflict.',u.eastFactory);
reject('AL04','Ammonium lactate','Price denominator not explicit; gross package weight 2 kg conflicts with 25 / 200 kg drum options.',u.kejian);
reject('AC01 / AC05 / AC09–AC11','Ammonium citrate','Triammonium citrate products excluded; scope is diammonium citrate only.',u.randa);
reject('AC02–AC04','Ammonium citrate','Conflicting CAS numbers and triammonium identity.',u.kolod);
reject('AC07','Diammonium citrate','Suspicious US$1/bag with no bag mass and multiple placeholder attributes.',u.essence);
reject('AC08','Diammonium citrate','Price header has no denominator; category /Ton and sample /kg do not establish a reliable bulk quote.',u.eastCitrate);
reject('New candidate: Mikem','Maleic acid','No numerical purity; product description contains an unrelated street-light template.','https://mikemchem.en.made-in-china.com/product/PpRYiAydaVWE/China-Top-Quality-Maleic-Acid-with-CAS-110-16-7-for-Organic-Intermediate.html');
reject('New candidate: Jinyu','Ammonium lactate','Attributes CAS 515-98-0 conflicts with description CAS 515-03-7.','https://qdjychemical.en.made-in-china.com/product/wOStPxUYHeaH/China-Factory-Supply-Ammonium-Lactate-Solution-CAS-515-98-0.html');
reject('New candidate: Ebos','Ammonium lactate','99% purity shown, but displayed quote has no explicit price denominator; not inferred from MOQ.','https://www.alibaba.com/product-detail/Amonium-Ammonium-Lactate-Cas-515-98_1600623185789.html');
reject('New candidate: Natural Micron','Ammonium lactate','60% concentration shown, but quote denominator and net package mass are not explicit.','https://www.alibaba.com/product-detail/Daily-Chemicals-Liquid-Ammonium-Lactate-CAS_1601587829339.html');
reject('New candidate: Leader','Ammonium lactate','Copied / implausible specifications: CAS in molecular-formula field, melting ≥350 °C, dark-green solubility test; insufficient product-specific evidence.','https://leaderbiogroup.lookchem.com/products/CasNo-515-98-0-China-Biggest-factory-Supply-High-Quality-AMMONIUM-LACTATE-CAS-515-98-0-27833383.html');
reject('New candidate: Wenxi','Ammonium lactate','Unrelated loratadine photographs / text and CAS in molecular-formula field; copied specifications excluded.','https://zzwenxi.lookchem.com/products/CasNo-515-98-0-Manufacturer-supply-CAS-515-98-0-with-best-quality-25826508.html');
reject('New candidate: Honghui','Ammonium lactate','60% bulk solution, but no numeric USD quotation.','https://www.ec21.com/product-details/Ammonium-Lactate-Solution-60-Food--9891396.html');
reject('New candidate: East Chemsources','Diammonium citrate','USD10 per metric ton appears to be a placeholder; excluded.','https://echemsouces.lookchem.com/products/CasNo-3012-65-5-AMMONIUM-CITRATE-diammoniumcitrate--30253175.html');
reject('New candidate: Kaiyu','Diammonium citrate','CAS and formula name dibasic citrate, but molecular weight 243.2 belongs to a different ammonium stoichiometry; specification identity uncertain.','https://www.chemicalbook.com/ProductDetail_EN_diammonium-citrate_3496191.htm');
reject('New candidate: Kolod','Diammonium citrate','Separate dibasic page has consistent CAS and packaging, but repeated $25 quotes across unrelated salts suggest template pricing; currency not explicitly USD.','https://www.chemicalbook.com/ProductDetail_EN_420641.htm');
reject('New candidate: Chuanghai','Diammonium citrate','Package and purity explicit, but dollar sign alone does not explicitly confirm USD; not retained under USD-only screening.','https://www.chemicalbook.com/ProductDetail_EN_1490063.htm');
reject('Search renewed: monoammonium','Ammonium maleate','CAS 44742-89-4. Directory displays prices, but no matching supplier product page with independently verifiable numeric USD quote found. Directory prices excluded.','https://www.guidechem.com/encyclopedia/ammonium-hydrogen-maleate-dic116530.html');
reject('Equation offer','Ammonium maleate','CAS 44742-89-4, purity and bulk drums shown; primary trade page requests latest price without numeric quote.','https://www.guidechem.com/trade/ammonium-hydrogen-maleate-id7581688.html');
reject('Search renewed: diammonium','Ammonium maleate','CAS 23705-99-9. Product page states 96% purity but has no public numeric quote; no qualified bulk price recorded.','https://www.alfa-chemistry.com/diammonium-maleate-cas-23705-99-9-item-108277.htm');

reject('Additional search: Pengjing','Ammonium lactate','20% / 60% solution concentrations explicit, but primary price header has no denominator and package mass is unspecified; excluded.','https://www.alibaba.com/product-detail/Professional-Supplier-Ammonium-Lactate-Solution-515_1601596747958.html');
reject('East Chemsources lactate','Ammonium lactate','Header quotes USD18–20/Kilogram but quantity tiers quote the same amounts per Metric Ton. Unit conflict and implausible total nitrogen 0.19% for claimed 99% product; excluded.','https://echemsouces.lookchem.com/products/CasNo-515-98-0-Ammonium-Lactate-30760637.html');
reject('Additional search: Top Pharm','Ammonium lactate','Page has 99% purity and 25 kg drums, but broad USD1–10/kg quote without quantity breaks appears generic; withheld under placeholder-price screening.','https://top-pharmchem.lookchem.com/products/CasNo-515-98-0-Ammonium-lactate-1585634.html');
reject('Additional search: Wuhan Pharmchem','Diammonium citrate','Supplier product link returns Product Not Available; no current quote verified.','https://www.alibaba.com/product-detail/Diammonium-strong-citrate-strong-Ammonium-strong_60650024845.html');

const wb=Workbook.create();
const review=wb.worksheets.add('Review'),listings=wb.worksheets.add('Listings'),props=wb.worksheets.add('Properties'),search=wb.worksheets.add('Search log');
const heads=['Record ID','Chemical','Chemical form','CAS as listed','Purity / concentration as listed','Quoted price (USD)','Currency','Quoted price unit','Overall MOQ value','Overall MOQ unit','Quantity required for quoted price','Packaging','Supplier','Verification method','Collected on','Selection','Review flag','Issues / notes','Supplier product URL','Unit-supporting URL'];
listings.getRange('A1').write([heads,...rows]);listings.tables.add(`A1:T${rows.length+1}`,true,'PriceListings');
listings.getRange(`F2:F${rows.length+1}`).setNumberFormat('0.00');listings.getRange(`O2:O${rows.length+1}`).setNumberFormat('yyyy-mm-dd');
listings.getRange(`P2:P${rows.length+1}`).dataValidation={rule:{type:'list',values:['Pending','Include','Exclude']}};listings.getRange(`P2:P${rows.length+1}`).format.fill='#FFF2CC';
props.getRange('A1').write([['Record ID','Property','Low value','High value','Unit','Conditions / chemical basis','Evidence / limitation','Source URL','Collected on'],...properties]);props.tables.add(`A1:I${properties.length+1}`,true,'ProductProperties');props.getRange(`I2:I${properties.length+1}`).setNumberFormat('yyyy-mm-dd');
search.getRange('A1').write([['Previous ID / candidate','Chemical','Disposition','Reason / search finding','Source URL','Reviewed on'],...log]);search.tables.add(`A1:F${log.length+1}`,true,'SearchLog');search.getRange(`F2:F${log.length+1}`).setNumberFormat('yyyy-mm-dd');
review.getRange('A2').values=[['Bulk chemical price collection']];review.getRange('A4:B4').values=[['Collected / reviewed on',date]];review.getRange('B4').setNumberFormat('yyyy-mm-dd');
review.getRange('A6:C10').values=[['Chemical','Product listings','Coverage note'],['Maleic acid',null,'6 suppliers; adds Surest and Sinobio; Xima tiers consolidated'],['Ammonium lactate',null,'4 suppliers; adds Youze and Sankai; expired status and price outlier flagged'],['Diammonium citrate',null,'4 suppliers; adds Dingyan, Dayang and GIHI; AC12 remains indexed-only'],['Ammonium maleate',null,'0 qualified quotes after renewed mono- and diammonium searches']];
for(let r=7;r<=10;r++)review.getRange(`B${r}`).formulas=[[`=COUNTIF(Listings!$B$2:$B$${rows.length+1},A${r})`]];
const notes=[
 ['Collection scope','Notes'],
 ['Price selection','Only lower bounds and lowest published quantity tiers are retained. An advertised low price is a screening proxy, not evidence of a negotiated bulk discount.'],
 ['USD / price units','Explicit USD and price denominator required. MOQ never supplies a missing price unit. Some Made-in-China denominators are corroborated by exact matching marketplace offers.'],
 ['Purity / concentration','Requires numerical purity, assay or concentration specific to the chemical. Generic active-substance content is excluded. Dry basis / solution basis remain as published, with gaps flagged.'],
 ['Identity / packaging','Conflicting CAS, conflicting package claims, triammonium citrate, apparent placeholder prices and copied specifications are excluded. Alternative package sizes are not themselves conflicts.'],
 ['Verification limits','13 product pages opened directly. AC12 uses indexed supplier-page text because live access presents bot verification. AL10 supplier page is marked expired; AL11 is a high price outlier. Displayed quotes do not certify current availability.'],
 ['Selection','All rows remain Pending. Confirm purity basis and the order quantity required for the quoted price before selecting inputs. Small supplier MOQs do not establish bulk discounts.'],
 ['Properties','Seller-reported properties remain tied to retained IDs. Unknown density units are not inferred. Density discrepancies may be retained with notes; properties from rejected products are removed.'],
 ['No conversions yet','No USD/kg conversion, dry-acid equivalent price, pooled market range or MSP calculation has been performed. Original quoted units and packaging are separate.'],
 ['Coverage gaps','No qualifying ammonium maleate quote was found. Missing fields are not filled from unrelated listings. Search log records exclusions and quote-only leads.'],
 ['Collection date','2026-10-31 is the access/review date. Indexed content can predate collection; AC12 was indexed approximately two months earlier.']
];
review.getRange('A13').write(notes);
for(const sh of [review,listings,props,search]){sh.showGridLines=false;sh.getUsedRange().format.font.name='Arial';sh.getUsedRange().format.font.size=10;sh.getUsedRange().format.verticalAlignment='center';}
function header(sh,r){const x=sh.getRange(r);x.format.fill='#284B63';x.format.font.color='#FFFFFF';x.format.font.bold=true;x.format.wrapText=true;x.format.rowHeight=40;}
review.getRange('A2').format.font.size=15;review.getRange('A2').format.font.bold=true;
review.getRange('A:A').format.columnWidth=29;review.getRange('B:B').format.columnWidth=19;review.getRange('C:C').format.columnWidth=80;
// Notes occupy a wide merged text area, while the flat data tables remain unmerged.
for(let r=13;r<=23;r++){review.getRange(`B${r}:C${r}`).merge();review.getRange(`B${r}:C${r}`).format.wrapText=true;}
review.getRange('A14:C23').format.rowHeight=58;review.getRange('A7:C10').format.rowHeight=42;review.getRange('C7:C10').format.wrapText=true;header(review,'A6:C6');header(review,'A13:C13');
for(const sh of [listings,props,search]){sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(2);sh.getUsedRange().format.rowHeight=86;sh.getUsedRange().format.wrapText=true;}
listings.getRange('A:A').format.columnWidth=12;listings.getRange('B:D').format.columnWidth=24;listings.getRange('E:E').format.columnWidth=43;listings.getRange('F:J').format.columnWidth=19;listings.getRange('K:M').format.columnWidth=38;listings.getRange('N:N').format.columnWidth=44;listings.getRange('O:P').format.columnWidth=17;listings.getRange('Q:Q').format.columnWidth=34;listings.getRange('R:R').format.columnWidth=72;listings.getRange('S:T').format.columnWidth=65;
props.getRange('A:A').format.columnWidth=12;props.getRange('B:B').format.columnWidth=24;props.getRange('C:E').format.columnWidth=18;props.getRange('F:G').format.columnWidth=52;props.getRange('H:H').format.columnWidth=65;props.getRange('I:I').format.columnWidth=16;
search.getRange('A:A').format.columnWidth=30;search.getRange('B:C').format.columnWidth=24;search.getRange('D:D').format.columnWidth=88;search.getRange('E:E').format.columnWidth=65;search.getRange('F:F').format.columnWidth=16;
header(listings,'A1:T1');header(props,'A1:I1');header(search,'A1:F1');
wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'Review!A6:C10',include:'values,formulas',tableMaxRows:5,tableMaxCols:3})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#N/A|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20}})).ndjson);
const expectedCAS={'Maleic acid':'110-16-7','Ammonium lactate':'515-98-0','Diammonium citrate':'3012-65-5'};
if(rows.some(r=>r[6]!=='USD'||!['kg','Kilogram','Ton'].includes(r[7])||!r[4].includes('%')||/active substances/i.test(r[4])||r[3]!==expectedCAS[r[1]]||r[5]<=0||!r[18].startsWith('https://'))||new Set(rows.map(r=>r[0])).size!==rows.length)throw Error('Listing validation failed');
if(properties.some(p=>!rows.some(r=>r[0]===p[0])))throw Error('Orphan property');
await fs.mkdir(previews,{recursive:true});
for(const [sheet,range,name] of [['Review','A1:C23','review'],['Listings','A10:F15','listings'],['Listings','F10:L15','prices'],['Listings','M10:R15','supplier-notes'],['Properties','A11:G17','properties'],['Search log','A24:D27','search-log']]){const img=await wb.render({sheetName:sheet,range,scale:1.4});await fs.writeFile(path.join(previews,`${name}.png`),new Uint8Array(await img.arrayBuffer()));}
const target=path.join(out,'bulk_chemical_market_prices.xlsx');await(await SpreadsheetFile.exportXlsx(wb)).save(target);await fs.copyFile(target,path.join(previews,'bulk_chemical_market_prices.xlsx'));
console.log(JSON.stringify({records:rows.length,properties:properties.length,searchLog:log.length,path:target}));
