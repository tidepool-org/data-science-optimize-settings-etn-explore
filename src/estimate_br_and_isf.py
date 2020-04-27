"""
Author: Ed Nykaza
Created: 2020-03-13
Last Updated: 2020-04-26

This is the first of several metrics that I will be building to help with determining Loop Settings.
Please note:
*   This code has not been peer-reviewed yet and should not be used to make medical decisions.
*   This settings estimator assumes that the user is using DIY Loop.

Here is the general model that is being solved using bounded or constrained optimization:

dBG[t,v] = EGP[t] - insulinEffect[t] + carbEffect[t] + activityEffect[t] + timingEffects[t] + valueEffects[t,v]

Though, for this version of the code, only the first two factors are being solved for:
*   dBG[t,v] = EGP[t] - insulinEffect[t]
*   dBG[t,v] = ISF(basalRate/60 - insulinEffect[t])

----------
CHANGE LOG:
    * only require Tidepool username (email) and password


"""

# %% REQUIRED LIBRARIES
import os
import sys
import getpass
import numpy as np
import pandas as pd
import datetime as dt
from pytz import timezone
from scipy.optimize import curve_fit, brute, fmin
import plotly.graph_objs as go
from plotly.offline import plot
import plotly.express as px
import subprocess
from dotenv import load_dotenv, find_dotenv

# %% TIDEPOOL API
if not os.path.exists("donor-data-pipeline"):
    process = subprocess.Popen(
        [
            "git",
            "clone",
            "-b",
            "jam/tidepool-api-improvements",
            "--single-branch",
            "https://github.com/tidepool-org/donor-data-pipeline.git",
            "donor-data-pipeline",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    output, errors = process.communicate()
    output = output.decode("utf-8")
    errors = errors.decode("utf-8")

sys.path.append(os.path.join("donor-data-pipeline", "src"))
import get_single_tidepool_dataset  # noqa: E402


# %% LOAD IN LOCAL ENV FILE (IF IT EXISTS)
# find .env automatically by walking up directories until it's found
dotenv_path = find_dotenv()

# load up the entries as environment variables
load_dotenv(dotenv_path)

# %% CONSTANTS
EPS = sys.float_info.epsilon
MGDL_PER_MMOLL = 18.01559


# %% GLOBAL FUNCTIONS
def getStartAndEndTimes(df, dateTimeField):
    dfBeginDate = df[dateTimeField].min()
    dfEndDate = df[dateTimeField].max()

    return dfBeginDate, dfEndDate


def removeDuplicates(df, criteriaDF):
    nBefore = len(df)
    df = df.loc[~(df[criteriaDF].duplicated())]
    df = df.reset_index(drop=True)
    nDuplicatesRemoved = nBefore - len(df)

    return df, nDuplicatesRemoved


def round_time(
    df,
    timeIntervalMinutes=5,
    timeField="time",
    roundedTimeFieldName="roundedTime",
    startWithFirstRecord=True,
    verbose=False,
):
    """
    A general purpose round time function that rounds the "time"
    field to nearest <timeIntervalMinutes> minutes
    INPUTS:
        * a dataframe (df) that contains a time field that you want to round
        * timeIntervalMinutes (defaults to 5 minutes given that most cgms output every 5 minutes)
        * timeField to round (defaults to the UTC time "time" field)
        * roundedTimeFieldName is a user specified column name (defaults to roundedTime)
        * startWithFirstRecord starts the rounding with the first record if True,
        and the last record if False (defaults to True)
        * verbose specifies whether the extra columns used to make calculations are returned
    """

    df.sort_values(by=timeField, ascending=startWithFirstRecord, inplace=True)
    df.reset_index(drop=True, inplace=True)

    # make sure the time field is in the right form
    t = pd.to_datetime(df[timeField])

    # calculate the time between consecutive records
    t_shift = pd.to_datetime(df[timeField].shift(1))
    df["timeBetweenRecords"] = (
        round(
            (t - t_shift).dt.days * (86400 / (60 * timeIntervalMinutes))
            + (t - t_shift).dt.seconds / (60 * timeIntervalMinutes)
        )
        * timeIntervalMinutes
    )

    # separate the data into chunks if timeBetweenRecords is greater than
    # 2 times the <timeIntervalMinutes> minutes so the rounding process starts over
    largeGaps = list(df.query("abs(timeBetweenRecords) > " + str(timeIntervalMinutes * 2)).index)
    largeGaps.insert(0, 0)
    largeGaps.append(len(df))

    for gIndex in range(0, len(largeGaps) - 1):
        chunk = t[largeGaps[gIndex] : largeGaps[gIndex + 1]]
        firstRecordChunk = t[largeGaps[gIndex]]

        # calculate the time difference between each time record and the first record
        df.loc[largeGaps[gIndex] : largeGaps[gIndex + 1], "minutesFromFirstRecord"] = (
            chunk - firstRecordChunk
        ).dt.days * (86400 / (60)) + (chunk - firstRecordChunk).dt.seconds / (60)

        # then round to the nearest X Minutes
        # NOTE: the ".000001" ensures that mulitples of 2:30 always rounds up.
        df.loc[largeGaps[gIndex] : largeGaps[gIndex + 1], "roundedMinutesFromFirstRecord"] = round(
            (df.loc[largeGaps[gIndex] : largeGaps[gIndex + 1], "minutesFromFirstRecord"] / timeIntervalMinutes)
            + 0.000001
        ) * (timeIntervalMinutes)

        roundedFirstRecord = (firstRecordChunk + pd.Timedelta("1microseconds")).round(str(timeIntervalMinutes) + "min")
        df.loc[largeGaps[gIndex] : largeGaps[gIndex + 1], roundedTimeFieldName] = roundedFirstRecord + pd.to_timedelta(
            df.loc[largeGaps[gIndex] : largeGaps[gIndex + 1], "roundedMinutesFromFirstRecord"], unit="m",
        )

    # sort by time and drop fieldsfields
    df.sort_values(by=timeField, ascending=startWithFirstRecord, inplace=True)
    df.reset_index(drop=True, inplace=True)
    if verbose is False:
        df.drop(
            columns=["timeBetweenRecords", "minutesFromFirstRecord", "roundedMinutesFromFirstRecord"], inplace=True,
        )

    return df


def convert_to_local_time(utc_time, current_timezone):
    if pd.notnull(utc_time):
        local_dt = utc_time.tz_convert(timezone(current_timezone))
    else:
        local_dt = np.nan
    return local_dt


def remove_timezone(local_timezoneAware):
    local_time = local_timezoneAware.tz_localize(None)
    return local_time


def add_uploadDateTime(df):
    if "upload" in data.type.unique():
        uploadTimes = pd.DataFrame(df[df.type == "upload"].groupby("uploadId").time.describe()["top"])
    else:
        uploadTimes = pd.DataFrame(columns=["top"])
    # if an upload does not have an upload date, then add one
    # NOTE: this is a new fix introduced with healthkit data...we now have
    # data that does not have an upload record
    unique_uploadIds = set(df["uploadId"].unique())
    unique_uploadRecords = set(df.loc[df["type"] == "upload", "uploadId"].unique())
    uploadIds_missing_uploadRecords = unique_uploadIds - unique_uploadRecords

    for upId in uploadIds_missing_uploadRecords:
        last_upload_time = df.loc[df["uploadId"] == upId, "time"].max()
        uploadTimes.loc[upId, "top"] = last_upload_time

    uploadTimes.reset_index(inplace=True)
    uploadTimes.rename(columns={"top": "uploadTime", "index": "uploadId"}, inplace=True)
    df = pd.merge(df, uploadTimes, how="left", on="uploadId")
    df["uploadTime"] = pd.to_datetime(df["uploadTime"])

    return df


def tempRemoveFields(df, removeFields):
    tempRemoveFields = list(set(df) & set(removeFields))
    tempDf = df[tempRemoveFields]
    df = df.drop(columns=tempRemoveFields)

    return df, tempDf


def flattenJson(df, doNotFlattenList):
    # remove fields that we don't want to flatten
    df, holdData = tempRemoveFields(df, doNotFlattenList)

    # get a list of data types of column headings
    columnHeadings = list(df)

    # loop through each columnHeading
    newDataFrame = pd.DataFrame()

    for colHead in columnHeadings:
        if any(isinstance(item, list) for item in df[colHead]):
            listBlob = df[colHead][df[colHead].astype(str).str[0] == "["]
            df.loc[listBlob.index, colHead] = df.loc[listBlob.index, colHead].str[0]

        # if the df field has embedded json
        if any(isinstance(item, dict) for item in df[colHead]):
            # grab the data that is in brackets
            jsonBlob = df[colHead][df[colHead].astype(str).str[0] == "{"]

            # replace those values with nan
            df.loc[jsonBlob.index, colHead] = np.nan

            # turn jsonBlob to dataframe
            newDataFrame = pd.concat(
                [newDataFrame, pd.DataFrame(jsonBlob.tolist(), index=jsonBlob.index).add_prefix(colHead + ".")], axis=1,
            )

    df = pd.concat([df, newDataFrame, holdData], axis=1)
    df.sort_index(axis=1, inplace=True)

    return df


def removeNegativeDurations(df):
    if "duration" in list(df):

        nNegativeDurations = sum(df.duration < 0)
        if nNegativeDurations > 0:
            df = df[~(df.duration < 0)]
    else:
        nNegativeDurations = np.nan

    return df, nNegativeDurations


def removeInvalidCgmValues(df):
    nBefore = len(df)
    # remove values < 38 and > 402 mg/dL
    df = df.drop(df[((df.type == "cbg") & (df.value < 2.109284236597303))].index)
    df = df.drop(df[((df.type == "cbg") & (df.value > 22.314006924003046))].index)
    nRemoved = nBefore - len(df)

    return df, nRemoved


def tslimCalibrationFix(df):
    if "payload.calibration_reading" in list(df):

        searchfor = ["tan"]
        tandemDataIndex = (df.deviceId.str.contains("|".join(searchfor))) & (df.type == "deviceEvent")

        payloadCalReadingIndex = df["payload.calibration_reading"].notnull()

        nTandemAndPayloadCalReadings = sum(tandemDataIndex & payloadCalReadingIndex)

        if nTandemAndPayloadCalReadings > 0:
            # if reading is > 30 then it is in the wrong units
            if df["payload.calibration_reading"].min() > 30:
                df.loc[payloadCalReadingIndex, "value"] = (
                    df[tandemDataIndex & payloadCalReadingIndex]["payload.calibration_reading"] / MGDL_PER_MMOLL
                )
            else:
                df.loc[payloadCalReadingIndex, "value"] = df[tandemDataIndex & payloadCalReadingIndex][
                    "payload.calibration_reading"
                ]
    else:
        nTandemAndPayloadCalReadings = 0

    return df, nTandemAndPayloadCalReadings


def removeCgmDuplicates(df, timeCriterion):
    if timeCriterion in df:
        df.sort_values(
            by=[timeCriterion, "uploadTime"], ascending=[False, False], inplace=True,
        )
        dfIsNull = df[df[timeCriterion].isnull()]
        dfNotNull = df[df[timeCriterion].notnull()]
        dfNotNull, nDuplicatesRemoved = removeDuplicates(dfNotNull, [timeCriterion, "value"])
        df = pd.concat([dfIsNull, dfNotNull])
        df.sort_values(
            by=[timeCriterion, "uploadTime"], ascending=[False, False], inplace=True,
        )
    else:
        nDuplicatesRemoved = 0

    return df, nDuplicatesRemoved


def mmolL_to_mgdL(mmolL):
    return mmolL * MGDL_PER_MMOLL


def clean_data(df):
    metadata = pd.DataFrame(index=["cleaning_data"])
    df["utcTimezoneAware"] = pd.to_datetime(df["time"], utc=True)
    df["utcTime"] = df["utcTimezoneAware"].dt.tz_convert(None)

    # get timezone from healthkit df (if it exists)
    for i in df.loc[df.payload.notnull(), "payload"].keys():
        if "HKTimeZone" in df["payload"][i].keys():
            df.loc[i, "timezone"] = df["payload"][i]["HKTimeZone"]

    # quick fix if you are missing timezone df
    if "timezone" not in list(df):
        manual_timezone = input(
            "Your timezone is missing (possible backend issue?),\n"
            + "please enter in the timezone that you spent the last 4 weeks in\n"
            + "(e.g., America/Chicago with no quotes, the format has to be perfect):\n"
            + "https://en.wikipedia.org/wiki/List_of_tz_database_time_zones \n"
        )
        df["timezone"] = manual_timezone

    # estimate local time (simple method)
    # TODO: use the latest estiamte local time code
    df.sort_values(by="time", ascending=False, inplace=True)
    df["timezone"].fillna(method="ffill", inplace=True)
    df["timezone"].fillna(method="bfill", inplace=True)
    df["localTimezoneAware"] = df[["utcTimezoneAware", "timezone"]].apply(lambda x: convert_to_local_time(*x), axis=1)
    df["localTime"] = df["localTimezoneAware"].apply(remove_timezone)

    # ADD UPLOAD DATE
    data_w_uploadTime = add_uploadDateTime(df)

    # FLATTEN COLUMNS THAT HAVE EMBEDDED JSON
    do_not_flatten_list = ["suppressed", "recommended", "payload"]
    flat_data = flattenJson(data_w_uploadTime, do_not_flatten_list)

    # CLEAN df
    # remove negative durations
    flat_data["duration"] = flat_data["duration"].astype(float)
    clean_data, nNegativeDurations = removeNegativeDurations(flat_data)
    metadata["nNegativeDurations"] = nNegativeDurations

    # get rid of cgm values too low/high (< 38 & > 402 mg/dL)
    clean_data, nInvalidCgmValues = removeInvalidCgmValues(clean_data)
    metadata["nInvalidCgmValues"] = nInvalidCgmValues

    # Tslim calibration bug fix
    clean_data, nTandemAndPayloadCalReadings = tslimCalibrationFix(clean_data)
    metadata["nTandemAndPayloadCalReadings"] = nTandemAndPayloadCalReadings

    # round all df to the nearest 1 minutes
    clean_data = round_time(
        clean_data,
        timeIntervalMinutes=1,
        timeField="localTime",
        roundedTimeFieldName="roundedLocalTime",
        startWithFirstRecord=True,
        verbose=False,
    )

    print("after cleaning there are", len(clean_data), "rows of data")

    # GET CGM DATA
    # group data by type
    groupedData = clean_data.groupby(by="type")

    # filter by cgm and sort by uploadTime
    cgmData = groupedData.get_group("cbg").dropna(axis=1, how="all")

    # get rid of duplicates that have the same ["deviceTime", "value"]
    cgmData, nCgmDuplicatesRemovedDeviceTime = removeCgmDuplicates(cgmData, "deviceTime")
    metadata["nCgmDuplicatesRemovedDeviceTime"] = nCgmDuplicatesRemovedDeviceTime

    # get rid of duplicates that have the same ["time", "value"]
    cgmData, nCgmDuplicatesRemovedUtcTime = removeCgmDuplicates(cgmData, "time")
    metadata["cnCgmDuplicatesRemovedUtcTime"] = nCgmDuplicatesRemovedUtcTime

    # get rid of duplicates that have the same "roundedTime"
    cgmData, nCgmDuplicatesRemovedRoundedTime = removeDuplicates(cgmData, "roundedLocalTime")
    metadata["nCgmDuplicatesRemovedRoundedTime"] = nCgmDuplicatesRemovedRoundedTime

    # get start and end times
    cgmBeginDate, cgmEndDate = getStartAndEndTimes(cgmData, "roundedLocalTime")
    metadata["cgm.beginDate"] = cgmBeginDate
    metadata["cgm.endDate"] = cgmEndDate

    # create a contiguous time series
    rng = pd.date_range(cgmBeginDate, cgmEndDate, freq="1min")
    contiguousData = pd.DataFrame(rng, columns=["cDateTime"])

    # get data in mg/dL units
    cgmData["mg_dL"] = mmolL_to_mgdL(cgmData["value"]).astype(int)

    print("rounding data to the nearest 1 minute and creating a 1 minute contiguous time series")

    # merge data
    contig_cgm = pd.merge(
        contiguousData,
        cgmData[["roundedLocalTime", "mg_dL"]],
        left_on="cDateTime",
        right_on="roundedLocalTime",
        how="left",
    )
    print("filling gaps in the cgm data and smoothing to 1-minute resolution")
    # fill gaps in cgm data if there are at most 2 hours of missing records in a row
    contig_cgm["mg_dL"].interpolate(limit=120, inplace=True)

    # smooth cgm data

    # Use centered rolling average over 15 minute time period,  ± 7 points
    print("smoothing data using a 120 minute centered window")
    contig_cgm["mg_dL_smooth"] = contig_cgm["mg_dL"].rolling(window=120, min_periods=1, center=True).mean()

    # sort descendingly
    contig_cgm = contig_cgm.sort_values(by="cDateTime", ascending=False).reset_index(drop=True)

    # calculate delta bg or bg rate
    contig_cgm["deltaBg"] = contig_cgm["mg_dL_smooth"] - contig_cgm["mg_dL_smooth"].shift(-1)

    # Use centered rolling average over 15 minute time period,  ± 7 points
    contig_cgm["deltaBg_smooth"] = contig_cgm["deltaBg"].rolling(window=120, min_periods=1, center=True).mean()

    contig_cgm["accelBg"] = contig_cgm["deltaBg_smooth"] - contig_cgm["deltaBg_smooth"].shift(-1)

    # Use centered rolling average over 15 minute time period,  ± 7 points
    contig_cgm["accelBg_smooth"] = contig_cgm["accelBg"].rolling(window=120, min_periods=1, center=True).mean()

    contig_cgm["jerkBg"] = contig_cgm["accelBg_smooth"] - contig_cgm["accelBg_smooth"].shift(-1)

    contig_cgm.drop(columns="roundedLocalTime", inplace=True)
    print("there are", len(contig_cgm), "cgm data points")

    # GET BASAL DATA
    # TODO: double check the frequency at which different insulin pumps
    # deliver basals and boluses to better capture the actual delivery amounts
    # per minute.
    # query data, drop null columns, and sort by time

    basal_data = clean_data[clean_data.type == "basal"].copy().dropna(axis=1, how="all")

    basal_data.sort_values("localTime", ascending=False, inplace=True)

    basal_data.reset_index(drop=True, inplace=True)
    basal_data[["type", "time", "localTime"]].head()

    # % PROCESS BASAL DATA

    # get rid of basal durations that are unrealistic
    nUnrealisticBasalDuration = (basal_data.duration < 0) | (basal_data.duration > 86400000)
    metadata["nUnrealisticBasalDuration"] = sum(nUnrealisticBasalDuration)
    basal_data.loc[nUnrealisticBasalDuration, "duration"] = np.nan

    # calculate the total amount of insulin delivered (duration * rate)
    basal_data.sort_values("localTime", ascending=True, inplace=True)
    basal_data.reset_index(inplace=True)

    basal_data["durationHours"] = basal_data["duration"] / 1000.0 / 3600.0
    basal_data["durationMinutes"] = basal_data["durationHours"] * 60
    basal_data["durationFromTimeMins"] = (
        (basal_data["localTime"].shift(-1) - basal_data["localTime"]).dt.seconds / 60
    ) + ((basal_data["localTime"].shift(-1) - basal_data["localTime"]).dt.days * 1440)
    basal_data["timeIssues"] = np.abs(basal_data["durationFromTimeMins"] - basal_data["durationMinutes"])
    basal_data["timeIssues"] = basal_data["timeIssues"].astype(float).round()
    basal_data["totalAmountOfBasalInsulin"] = basal_data["durationHours"] * basal_data["rate"]
    basal_data.head()

    # GET OTHER HEALTHKIT INFORMATION
    for i in basal_data.index:
        basal_data.loc[i, "basalDeliveryReason"] = basal_data.loc[i, "payload"]["HKInsulinDeliveryReason"]
        if pd.notnull(basal_data.loc[i, "suppressed"]):
            basal_data.loc[i, "scheduledBasalRate"] = basal_data.loc[i, "suppressed"]["rate"]

    print("there are", len(basal_data), "basal insulin events")

    # CREATE A CONTINGUOUS BASAL TIME SERIES
    basal_30 = basal_data.copy()

    # % ROUND BASALS TO THE NEAREST 30 seconds
    basal_30["ratePer30sec"] = basal_30["rate"] / 120

    basal_30 = round_time(
        basal_30,
        timeIntervalMinutes=0.5,
        timeField="localTime",
        roundedTimeFieldName="roundedLocalTime",
        startWithFirstRecord=True,
        verbose=False,
    )

    # % CREATE A CONTIGUOUS BASAL TIME SERIES
    basalBeginDate, basalEndDate = getStartAndEndTimes(basal_30, "roundedLocalTime")
    metadata["basal.beginDate"] = basalBeginDate
    metadata["basal.endDate"] = basalEndDate

    rng = pd.date_range(basalBeginDate, basalEndDate, freq="0.5min")
    contiguousData = pd.DataFrame(rng, columns=["cDateTime"])

    # merge data
    temp_contig_basal = pd.merge(
        contiguousData, basal_30, left_on="cDateTime", right_on="roundedLocalTime", how="left",
    )

    # fill forward in the missing basal rates
    temp_contig_basal["ratePer30sec"].fillna(method="ffill", inplace=True)

    # resample to get the rate per 1 minutes
    series = pd.Series(
        temp_contig_basal["ratePer30sec"].values, name="ratePer30sec", index=temp_contig_basal["cDateTime"],
    )

    contig_basal = pd.DataFrame(series.resample("1min").sum()).reset_index()
    contig_basal.rename(columns={"ratePer30sec": "basalRatePerMin"}, inplace=True)
    contig_basal["basalRatePerHour"] = contig_basal["basalRatePerMin"] * 60

    # add the scheduled basal rate and other healthkit data back in
    sbr = basal_data[["roundedLocalTime", "basalDeliveryReason", "scheduledBasalRate"]].copy()
    sbr, nDups = removeDuplicates(sbr, "roundedLocalTime")

    # process the bolus data
    bolus_data = clean_data[clean_data.type == "bolus"].dropna(axis=1, how="all")

    for i in bolus_data.index:
        bolus_data.loc[i, "bolusDeliveryReason"] = bolus_data.loc[i, "payload"]["HKInsulinDeliveryReason"]
    bolus_data[["roundedLocalTime", "normal", "bolusDeliveryReason"]]
    bolus_groups = bolus_data.groupby("roundedLocalTime")
    combined_boluses = pd.DataFrame(bolus_groups["normal"].sum()).reset_index()
    bolus_data, nDups = removeDuplicates(bolus_data, "roundedLocalTime")
    combined_boluses = pd.merge(
        combined_boluses, bolus_data[["roundedLocalTime", "bolusDeliveryReason"]], on="roundedLocalTime", how="left",
    )

    print("there are", len(bolus_data), "bolus insulin events")

    # merge the bolus records with the contiguous basal insulin records
    contig_insulin = pd.merge(
        contig_basal, combined_boluses, left_on="cDateTime", right_on="roundedLocalTime", how="left",
    )
    contig_insulin.drop(columns="roundedLocalTime", inplace=True)

    contig_insulin = pd.merge(contig_insulin, sbr, left_on="cDateTime", right_on="roundedLocalTime", how="left",)

    contig_insulin["scheduledBasalRate"].fillna(method="ffill", inplace=True)
    contig_insulin["scheduledBasalRate"].fillna(method="bfill", inplace=True)

    contig_insulin["insulinAmount"] = contig_insulin["basalRatePerMin"] + contig_insulin["normal"].fillna(0)
    contig_insulin.drop(columns="roundedLocalTime", inplace=True)
    print("combining basal and bolus records into an insulin time series")

    # GET CARB DATA

    carb_data = clean_data[clean_data.type == "food"].dropna(axis=1, how="all")
    for i in carb_data.index:
        carb_data.loc[i, "grams"] = carb_data.loc[i, "nutrition.carbohydrate"]["net"]
        for k in carb_data.loc[i, "payload"].keys():
            if "AbsorptionTimeMinutes" in k:
                absorb_key = k
                carb_data.loc[i, "carbAbsorbTimeMinutes"] = carb_data.loc[i, "payload"][absorb_key]

    # get rid of duplicates that have the same ["time", "value"]
    carb_data, ncarb_dataDuplicatesRemovedUtcTime = removeDuplicates(carb_data, "utcTime")
    metadata["ncarb_dataDuplicatesRemovedUtcTime"] = ncarb_dataDuplicatesRemovedUtcTime
    print("there are", len(carb_data), "carb events")

    # combine the all of the data together
    all_df = pd.merge(contig_insulin, contig_cgm, on="cDateTime", how="left")

    all_df = pd.merge(
        all_df,
        carb_data[["roundedLocalTime", "grams", "carbAbsorbTimeMinutes"]],
        left_on="cDateTime",
        right_on="roundedLocalTime",
        how="left",
    )
    all_df.drop(columns=["roundedLocalTime"], inplace=True)
    all_df["hour"] = all_df["cDateTime"].dt.hour
    all_df["day"] = all_df["cDateTime"].dt.date
    all_df["month"] = all_df["cDateTime"].dt.month
    all_df["year"] = all_df["cDateTime"].dt.year

    # infer the CIR given that it is not given in healthkit records
    all_df["cirInferred"] = np.round(all_df["grams"] / all_df["normal"], 1)

    # create a bolus amount for carbs, for plotting purposes
    median_cir = all_df.cirInferred.median()
    all_df["bolus_height"] = all_df["normal"]
    for i in all_df[(all_df["grams"].notnull()) & (all_df["normal"].isnull())].index:
        all_df.loc[i, "bolus_height"] = all_df.loc[i, "grams"] / median_cir

    all_df[(all_df["grams"].notnull()) & (all_df["normal"].isnull())]

    all_df.sort_values("cDateTime", inplace=True)
    all_df.reset_index(drop=True, inplace=True)

    print("combining all of the data into a single data frame")

    return all_df


def calc_loop_insulin_models(df):
    # calculate the expected drop in bg/time (bg rate) for the insulin model
    # start by usign the loop model
    t = np.arange(0, 6 * 60, 1)

    insulin_amount_matrix = get_insulin_amount_matrix(df["insulinAmount"].values, insulin_curve_len=len(t))

    # NOTE: a peak time of 75 mintues + 10 minute delay = 85 minute peak time
    for peak_time, model_name in zip([55, 65, 75], ["fiasp", "child", "adult"]):
        normalized_active_insulin = loop_exponential_model(t, peak_time_minutes=peak_time, delay=10)

        reshaped_normalized_active_insulin = np.reshape(normalized_active_insulin, (1, len(normalized_active_insulin)))

        dG_insulin_loop = np.dot(reshaped_normalized_active_insulin, insulin_amount_matrix)

        df["dG_insulin_loop_" + model_name] = dG_insulin_loop[0]

    return df


def constrained_optimization(
    function_to_fit, known_input_array, measured_output_array, initial_parameter_estimates, lower_bounds, upper_bounds,
):
    """
    Constrained optimaization using scipy curve_fit function.

    A wrapper or helper function to
    `scipy curve_fit <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_
    function given that the curve_fit function is a little esoteric.
    This function uses the trust region reflection algorithm to search for the
    optimal value of each parameter, and searchers over the parameter bounds.


    Parameters
    ----------
    function_to_fit : function
        The function you want to optimize over.
    known_input_array : numpy array
        The known inputs to the function.
    measured_output_array : numpy array
        The measured or estimated outputs to your function.
    initial_parameter_estimates : list(param1, param2, ..., paramN)
        The optimization will start with this value when searching
        for the optimal parameter value.
    lower_bounds : list(param1, param2, ..., paramN)
        The lower bound for each parameter.
    upper_bounds : list(param1, param2, ..., paramN)
        The upper bound for each parameter.

    Returns
    -------
    tuple
        model_goodness_of_fit_rmse : the closer to 0 the better
        optimal_parameter_values,
        parameter_standard_deviation,
        covariance_matrix

    Example
    -------
    Simple example that finds slope (m) and y-intercept (b) to a line.
    >>> def estimate_y(x, m, b):
            y_hat = m*x + b
            return y_hat

    >>> y = np.arange(0, 1000)
    >>> x = np.arange(0, len(y))
    >>> m_actual = 1
    >>> b_actual = 1
    >>> y_measured = (
            m_actual*x
            + b_actual
            + np.random.normal(loc=0, scale=0.0001, size=len(x))
        )

    >>> output_tuple = constrained_optimization(
            function_to_fit=estimate_y,
            known_input_array=x,
            measured_output_array=y_measured,
            initial_parameter_estimates=[0, 0],
            lower_bounds=[-np.inf, -np.inf],
            upper_bounds=[np.inf, np.inf],
        )
    >>> print("model_goodness_of_fit_rmse =", output_tuple[0])
    model_goodness_of_fit_rmse = 9.938141475168524e-05
    >>> print("optimal_parameter_values =", output_tuple[1])
    optimal_parameter_values = [1.         1.00000632]
    >>> print("parameter_standard_deviation =", output_tuple[2])
    parameter_standard_deviation = [1.08867106e-08 6.28073135e-06]
    >>> print("covariance_matrix =", output_tuple[3])
    covariance_matrix = [[ 1.18520467e-16 -5.92010068e-14]
     [-5.92010068e-14  3.94475863e-11]]
    """

    optimal_parameter_values, covariance_matrix = curve_fit(
        f=function_to_fit,
        xdata=known_input_array,
        ydata=measured_output_array,
        p0=initial_parameter_estimates,
        bounds=(lower_bounds, upper_bounds),
        check_finite=True,
        method="trf",
    )

    fit_residuals = measured_output_array - function_to_fit(known_input_array, *optimal_parameter_values)

    model_goodness_of_fit_rmse = (np.sum(fit_residuals ** 2) / (fit_residuals.size - 2)) ** 0.5

    parameter_standard_deviation = np.sqrt(np.diag(covariance_matrix))

    output = (
        model_goodness_of_fit_rmse,
        optimal_parameter_values,
        parameter_standard_deviation,
        covariance_matrix,
    )

    return output


def get_insulin_amount_matrix(insulin_amount_array, insulin_curve_len=96):
    T = len(insulin_amount_array)
    matrix_shape = (insulin_curve_len, T)
    insulin_amount_matrix = np.zeros(matrix_shape)
    rows, cols = np.indices(matrix_shape)
    for k in np.arange(0, T):
        row_vals = np.diag(rows, k=k)
        col_vals = np.diag(cols, k=k)
        insulin_amount_matrix[row_vals, col_vals] = insulin_amount_array[k]

    return insulin_amount_matrix


def brute_optimization(
    function_to_fit, *known_inputs, lower_bounds, upper_bounds, mesh_size, workers=1,
):
    """
    Brute force optimaization using scipy brute function.

    <more details here>


    Parameters
    ----------
    function_to_fit : function
        The function you want to optimize over.
    <fill in details>
    lower_bounds : list(param1, param2, ..., paramN)
        The lower bound for each parameter.
    upper_bounds : list(param1, param2, ..., paramN)
        The upper bound for each parameter.

    Returns
    -------
    tuple
        model_goodness_of_fit_rmse : the closer to 0 the better
        optimal_parameter_values,
        parameter_standard_deviation,
        covariance_matrix

    Example
    -------
    """

    rranges = []
    n_parameters = len(lower_bounds)
    for r in np.arange(0, n_parameters):
        step_size = (upper_bounds[r] - lower_bounds[r]) / (mesh_size - 1)
        rranges = np.append(rranges, slice(lower_bounds[r], upper_bounds[r] + step_size, step_size),)

    brute_results = brute(
        function_to_fit,
        rranges,
        args=known_inputs,
        full_output=True,
        finish=fmin,  # fmin will look for a local minimum
        workers=workers,
    )

    optimal_parameter_values = brute_results[0]
    model_goodness_of_fit_rmse = brute_results[1]
    iteration_parameters = brute_results[2]

    iteration_loss_scores = brute_results[3]
    iteration_results = pd.DataFrame(iteration_loss_scores.reshape([-1, 1]), columns=["loss"])

    if len(np.shape(iteration_parameters)) == 1:
        iteration_results["param{}".format(r)] = iteration_parameters
    else:
        for r in np.arange(0, n_parameters):
            iteration_results["param{}".format(r)] = iteration_parameters[r].reshape([-1, 1])

    output = (
        model_goodness_of_fit_rmse,
        optimal_parameter_values,
        iteration_results,
    )

    return output


# GENERIC CARB AND INSULIN MODELS IN CESCON (eq. 8.2)
def cescon(
    t, delay, tau,
):
    """
    generic insulin and carb model given in Cescon's Thesis (2013)

    This equation is taken from eqs. 8.3 and 8.4, which given the same model
    for the effect of insulin on glucose (8.3) and carbs (8.4). Given that
    these equations give the impact on glucose, the amount of difference
    between glucose amounts at time (t + delta_t) and t is taken to approximate
    the activity curves. The general equation is:
        G(s) = e^(-s*delay) * K / s(1 + s*tau)

    Taking the Laplace Transform gives:
        G(t) = K * (1 - e^(-(t-delay)/tau)) * Heaviside(t - delay)

    NOTE:
        * K is a gain factor (i.e., ISF for insulin and ISF/CIR for carbs),
        which is removed or set to 1 so that the gain can be applied elsewhere.
        * delay is written as tau in the Cescon reference
        * tau is written as T in the Cescon reference

    Parameters
    ----------
    t : array
        DESCRIPTION.
    delay : int
        DESCRIPTION.
    tau : float
        DESCRIPTION.

    Returns
    -------
    normalized_activity_curve : numpy array
        DESCRIPTION.

    """
    bg = (1 - (np.exp(-(t - delay) / tau))) * np.heaviside(t - delay, 1)
    activity_curve = np.append(0, np.diff(bg))
    normalized_activity_curve = activity_curve / np.sum(activity_curve)

    return normalized_activity_curve


def two_parameter_biexponential_activity_curve(
    t, tau1, tau2,  # 55,  # 70
):
    activity_curve = (1 / (tau2 - tau1)) * (np.exp(-t / tau2) - np.exp(-t / tau1))

    normalized_activity_curve = activity_curve / np.sum(activity_curve)

    return normalized_activity_curve


# loop insulin model
def loop_exponential_model(t, peak_time_minutes, delay):
    delay = int(delay)
    if delay < 2:
        print("delay must be >= 2, setting delay to 2")
        delay = 2

    pk = peak_time_minutes
    step_size = t[1] - t[0]
    dur = len(t) * step_size
    n_delay_zeros = np.int(delay / step_size)

    tau = pk * (1 - pk / dur) / (1 - 2 * pk / dur)
    a = 2 * tau / dur
    S = 1 / (1 - a + (1 + a) * np.exp(-dur / tau))

    iob_percent = 1 - S * (1 - a) * ((pow(t, 2) / (tau * dur * (1 - a)) - (t / tau) - 1) * np.exp(-t / tau) + 1)
    # this logic keeps peak at peak time + delay
    insulin_activity = np.diff(-iob_percent)
    activity_curve = np.append(np.zeros(n_delay_zeros), insulin_activity[: -(n_delay_zeros - 1)])
    normalized_activity_curve = activity_curve / np.sum(activity_curve)

    return normalized_activity_curve


# functions to get paratmers
def est_isf_br_and_egp_insulin_only_models(
    df,
    name="no_name_given",
    initial_parameter_estimates=[40, 1, 0.5],
    lower_bounds=[10, 0.025, 0.04],
    upper_bounds=[500, 4, 350],
):
    # TODO: make this more flexible by passing in models and using classes
    # TODO: if time step is NOT every minute, then basal rate
    fixed_egp = initial_parameter_estimates[2]

    def estimate_isf_fixed_egp(insulin_effect_array, isf, egp=fixed_egp):
        return egp - (isf * insulin_effect_array)

    def estimate_isf_and_br(insulin_effect_array, isf, br):
        return isf * ((br / 60) - (insulin_effect_array))

    def estimate_isf_and_egp(insulin_effect_array, isf, egp):
        return egp - (isf * insulin_effect_array)

    output_row_names = [
        "startDate",
        "endDate",
        "nMinutesOfData",
        "rmse",
        "isf",
        "isfSD",
        "br",
        "brSD",
        "egp",
        "egpSD",
    ]
    model_output = pd.DataFrame(index=output_row_names)
    start_date = df["cDateTime"].min()
    end_date = df["cDateTime"].max()
    N = len(df)
    model_output.loc[output_row_names, "all_models"] = (
        start_date,
        end_date,
        N,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    )

    for insulin_model in ["loop", "biexponential", "cescon"]:
        model_name = "est_isf_fixed_egp_{}".format(insulin_model)
        output_tuple = constrained_optimization(
            function_to_fit=estimate_isf_fixed_egp,
            known_input_array=df["dG_insulin_{}".format(insulin_model)].values,
            measured_output_array=df["deltaBg"].values,
            initial_parameter_estimates=initial_parameter_estimates[0],
            lower_bounds=lower_bounds[0],
            upper_bounds=upper_bounds[0],
        )
        # print(output_tuple)
        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)
        isf_sd = np.round(output_tuple[2][0], 2)
        br = np.round(fixed_egp * 60 / isf, 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            isf_sd,
            br,
            np.nan,
            fixed_egp,
            np.nan,
        )

        model_name = "est_isf_br_{}".format(insulin_model)
        output_tuple = constrained_optimization(
            function_to_fit=estimate_isf_and_br,
            known_input_array=df["dG_insulin_{}".format(insulin_model)].values,
            measured_output_array=df["deltaBg"].values,
            initial_parameter_estimates=initial_parameter_estimates[0:2],
            lower_bounds=lower_bounds[0:2],
            upper_bounds=upper_bounds[0:2],
        )

        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)
        isf_sd = np.round(output_tuple[2][0], 2)
        br = np.round(output_tuple[1][1], 2)
        br_sd = np.round(output_tuple[2][1], 2)
        egp = np.round(isf * br / 60, 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            isf_sd,
            br,
            br_sd,
            egp,
            np.nan,
        )

        model_name = "est_isf_egp_{}".format(insulin_model)
        output_tuple = constrained_optimization(
            function_to_fit=estimate_isf_and_egp,
            known_input_array=df["dG_insulin_{}".format(insulin_model)].values,
            measured_output_array=df["deltaBg"].values,
            initial_parameter_estimates=[initial_parameter_estimates[0], initial_parameter_estimates[2]],
            lower_bounds=[lower_bounds[0], lower_bounds[2]],
            upper_bounds=[upper_bounds[0], upper_bounds[2]],
        )

        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)
        isf_sd = np.round(output_tuple[2][0], 2)
        egp = np.round(output_tuple[1][1], 2)
        egp_sd = np.round(output_tuple[2][1], 2)
        br = np.round(egp * 60 / isf, 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            isf_sd,
            br,
            br_sd,
            egp,
            egp_sd,
        )

    output_df = pd.DataFrame(model_output.T.stack(), columns=[name]).T
    output_df.columns = [".".join(col).strip() for col in output_df.columns.values]

    return output_df


# RUN ON ALL DATA WITH BRUTE FORCE
def get_rmse(estimated, actual):
    return np.sqrt(np.mean((estimated - actual) ** 2))


def estimate_isf_fixed_egp_brute(params_to_estimate, knowns):
    isf = params_to_estimate
    egp, just_insulin_snippets, insulin_model = knowns
    insulin_effect = just_insulin_snippets["dG_insulin_{}".format(insulin_model)].values

    estimated = egp - (isf * insulin_effect)
    rmse_loss = get_rmse(estimated, just_insulin_snippets["deltaBg"].values)

    return rmse_loss


def estimate_isf_and_br_brute(params_to_estimate, knowns):
    isf, br = params_to_estimate
    just_insulin_snippets, insulin_model = knowns

    insulin_effect = just_insulin_snippets["dG_insulin_{}".format(insulin_model)].values

    estimated = isf * ((br / 60) - (insulin_effect))
    rmse_loss = get_rmse(estimated, just_insulin_snippets["deltaBg"].values)

    return rmse_loss


def estimate_isf_and_egp_brute(params_to_estimate, knowns):
    isf, egp = params_to_estimate
    just_insulin_snippets, insulin_model = knowns

    insulin_effect = just_insulin_snippets["dG_insulin_{}".format(insulin_model)].values

    estimated = egp - (isf * insulin_effect)
    rmse_loss = get_rmse(estimated, just_insulin_snippets["deltaBg"].values)

    return rmse_loss


def est_isf_br_and_egp_insulin_only_models_brute(
    df,
    name="brute_force_method",
    initial_parameter_estimates=[40, 1, 0.5],
    lower_bounds=[10, 0.025, 0.04],
    upper_bounds=[500, 4, 350],
):
    # TODO: make this more flexible by passing in models and using classes
    # TODO: if time step is NOT every minute, then basal rate
    fixed_egp = initial_parameter_estimates[2]

    output_row_names = [
        "startDate",
        "endDate",
        "nMinutesOfData",
        "rmse",
        "isf",
        "br",
        "egp",
    ]

    model_output = pd.DataFrame(index=output_row_names)
    start_date = df["cDateTime"].min()
    end_date = df["cDateTime"].max()
    N = len(df)
    model_output.loc[output_row_names, "all_models"] = (
        start_date,
        end_date,
        N,
        np.nan,
        np.nan,
        np.nan,
        np.nan,
    )

    for insulin_model in ["loop", "biexponential", "cescon"]:
        # find isf and with egp fixed
        model_name = "est_isf_fixed_egp_{}".format(insulin_model)
        known_inputs = (fixed_egp, df, insulin_model)
        output_tuple = brute_optimization(
            estimate_isf_fixed_egp_brute,
            known_inputs,
            lower_bounds=[lower_bounds[0]],
            upper_bounds=[upper_bounds[0]],
            mesh_size=5,
        )
        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            np.nan,
            fixed_egp,
        )

        # find isf and basal rate
        model_name = "est_isf_br_{}".format(insulin_model)
        known_inputs = (df, insulin_model)

        output_tuple = brute_optimization(
            estimate_isf_and_br_brute,
            known_inputs,
            lower_bounds=lower_bounds[0:2],
            upper_bounds=upper_bounds[0:2],
            mesh_size=3,
        )

        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)
        br = np.round(output_tuple[1][1], 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            br,
            np.nan,
        )

        # find isf and egp
        model_name = "est_isf_egp_{}".format(insulin_model)
        output_tuple = brute_optimization(
            estimate_isf_and_egp_brute,
            known_inputs,
            lower_bounds=[lower_bounds[0], lower_bounds[2]],
            upper_bounds=[upper_bounds[0], upper_bounds[2]],
            mesh_size=3,
        )
        model_rmse = np.round(output_tuple[0], 2)
        isf = np.round(output_tuple[1][0], 2)
        egp = np.round(output_tuple[1][1], 2)

        model_output.loc[output_row_names, model_name] = (
            start_date,
            end_date,
            N,
            model_rmse,
            isf,
            np.nan,
            egp,
        )

    output_df = pd.DataFrame(model_output.T.stack(), columns=[name]).T
    output_df.columns = [".".join(col).strip() for col in output_df.columns.values]

    return output_df


def get_isf_br_egp(brute_force_results):
    fixed_isf = np.round(
        np.average(
            [
                brute_force_results["est_isf_br_loop.isf"].values,
                brute_force_results["est_isf_br_biexponential.isf"].values,
                brute_force_results["est_isf_br_cescon.isf"].values,
                brute_force_results["est_isf_egp_loop.isf"].values,
                brute_force_results["est_isf_egp_biexponential.isf"].values,
                brute_force_results["est_isf_egp_cescon.isf"].values,
            ]
        ),
        0,
    )

    fixed_br = np.round(
        np.average(
            [
                brute_force_results["est_isf_br_loop.br"].values,
                brute_force_results["est_isf_br_biexponential.br"].values,
                brute_force_results["est_isf_br_cescon.br"].values,
            ]
        ),
        2,
    )

    fixed_egp = np.round(
        np.average(
            [
                brute_force_results["est_isf_egp_loop.egp"].values,
                brute_force_results["est_isf_egp_biexponential.egp"].values,
                brute_force_results["est_isf_egp_cescon.egp"].values,
            ]
        ),
        2,
    )

    return [fixed_isf, fixed_br, fixed_egp]


def estimate_insulin_curve_params_fixed_isf_br_brute(params_to_estimate, knowns):
    """parameters to estimate """
    param1, param2 = params_to_estimate
    """knowns, that are needed for the estimate """
    (insulin_model, t, insulin_amount_matrix, just_insulin_snippets, isf, br,) = knowns

    normalized_active_insulin = insulin_model(t, param1, param2)

    reshaped_normalized_active_insulin = np.reshape(normalized_active_insulin, (1, len(normalized_active_insulin)))

    insulin_effect = (np.dot(reshaped_normalized_active_insulin, insulin_amount_matrix))[0][
        just_insulin_snippets["index"].values
    ]

    estimated = isf * ((br / 60) - (insulin_effect))
    rmse_loss = get_rmse(estimated, just_insulin_snippets["deltaBg"].values)

    return rmse_loss


def prep_insulin_snippets(df):
    # only take cases where all of the cgm data and deltaBg is present
    df["hasCgm"] = df["deltaBg"].notnull()

    # 5 hours after the last carb
    just_carbs = df["grams"]
    carbs_within_last_5hours = just_carbs.rolling(60 * 5, min_periods=1).sum()
    df["carbsLast5Hour"] = carbs_within_last_5hours
    df["hasCarbsWithinLast5Hours"] = carbs_within_last_5hours > 0

    return df


def get_delta_bg_from_isf_br_loop_insulin_curve(insulin_effect_array, isf, br):
    return isf * ((br / 60) - (insulin_effect_array))


def find_just_insulin_snippets(
    df, field="justInsulin_adult", min_snippet_minutes=120,
):
    """
    Return array of snippet starts and sizes that meet input conditions.

    Parameters
    ----------
    df : pd.DataFrame
        Pandas dataframe that contains cgm data and a time field
    field : str
        Name of dataframe column that contains criterion
    min_snippet_minutes : int
        This defines the minimum length of a snippet in minutes. A cgm snippet must be
        greater than this threshold in order to be stitched together with another snippet.
        TODO: make this explanation clearer

    Returns
    -------
    np.array
        A two column array. The first column contains the index of the snippet starts,
        and the second column contains the size of each snippet.
    """

    # get the time when just insulin turns on (1) and off (-1)
    df["on_off"] = df[field] - df[field].shift(1)
    df["on_off"].fillna(0, inplace=True)

    # separate the data into chunks if timeBetweenRecords is greater than largest allowable gap
    snippet_start_index = list(df.query("on_off == 1").index)

    if len(snippet_start_index) > 0:
        temp_df = df[df["on_off"].index > snippet_start_index[0]]
        snippet_end_index = list(temp_df.query("on_off == -1").index - 1)
        if len(snippet_start_index) > len(snippet_end_index):
            snippet_end_index.append(df.index.max())

        # create a dataframe of just the snippet starts and sizes
        just_snippets = pd.DataFrame(snippet_start_index, columns=["snippet_start_index"])
        just_snippets["snippet_end_index"] = snippet_end_index
        just_snippets["snippet_size"] = just_snippets["snippet_end_index"] - just_snippets["snippet_start_index"]

        # get rid of all snippets less than (user defined length)
        keep_snippets_array = just_snippets.loc[just_snippets["snippet_size"] >= min_snippet_minutes, :].values

        just_insulin_snippets = pd.DataFrame(columns=df.columns)
        for i in range(0, len(keep_snippets_array)):
            t_df = df.loc[keep_snippets_array[i][0] : keep_snippets_array[i][1], :]
            just_insulin_snippets = pd.concat([just_insulin_snippets, t_df])

        just_insulin_snippets.reset_index(inplace=True)
        just_insulin_snippets["consec_start"] = just_insulin_snippets["index"].diff() != 1
        just_insulin_snippets["snippet_id"] = just_insulin_snippets["consec_start"].cumsum()
        total_snippets = just_insulin_snippets["consec_start"].sum()
    else:
        just_insulin_snippets = np.nan
        total_snippets = 0
    avg_snippets_per_week = total_snippets / weeks_of_data

    return just_insulin_snippets, avg_snippets_per_week


def get_insulin_snippets_given_isf_br(df, isf, br, snippet_threshold=1, insulin_model="adult"):
    df["dBg_non_insulin_" + insulin_model] = df["deltaBg_smooth"] - get_delta_bg_from_isf_br_loop_insulin_curve(
        df["dG_insulin_loop_" + insulin_model].values, isf, br
    )

    # find cases where the non-insulin BG rate is not increasing/decrease by more than +/- snippet_treshold mg/dL/min
    just_rates = df["dBg_non_insulin_" + insulin_model]
    rates_within_last_2hours = just_rates.abs().rolling(60 * 2, min_periods=1).max()
    df["hasRatesGtThresholdLast2Hours_" + insulin_model] = rates_within_last_2hours > snippet_threshold

    # get the individual time stamps that are just insulin and has cgm data
    df["justInsulin_" + insulin_model] = (
        ~(df["hasRatesGtThresholdLast2Hours_" + insulin_model] | df["hasCarbsWithinLast5Hours"])
        & (df["hasCgm"])
        & (df["hour"] < 12)  # add in constraint of the middle of the night
    )

    # get consecutive snippets where there is at least 2 hours of
    just_insulin_snippets, avg_snippets_per_week = find_just_insulin_snippets(df, field="justInsulin_" + insulin_model)

    return just_insulin_snippets, avg_snippets_per_week, df


def get_rmse_given_isf_br(isf, br, df, snippet_threshold, insulin_model):
    (insulin_snippets, avg_snippets_per_week, _,) = get_insulin_snippets_given_isf_br(
        df, isf, br, snippet_threshold=snippet_threshold, insulin_model=insulin_model,
    )

    if avg_snippets_per_week >= 3:
        delta_bg_est = get_delta_bg_from_isf_br_loop_insulin_curve(
            insulin_snippets["dG_insulin_loop_" + insulin_model].values, isf, br,
        )
        rmse = get_rmse(delta_bg_est, insulin_snippets["dG_insulin_loop_adult"].values)
        normalized_rmse = rmse / avg_snippets_per_week
    else:
        normalized_rmse = np.nan
    # print(insulin_model, isf, br, avg_snippets_per_week, normalized_rmse)
    return normalized_rmse


def make_evidence_plots(
    df, ji_df, snippet_threshold, insulin_model, days_to_show=7, title="None Given",
):
    END_DATE = df["cDateTime"].max()
    START_DATE = df["cDateTime"].max() - pd.Timedelta(days_to_show, "day")

    bg_smooth_trace = go.Scattergl(
        x=ji_df["cDateTime"],
        y=ji_df["mg_dL_smooth"],
        name="Just Insulin Snippets",
        mode="markers",
        opacity=1,
        marker=dict(size=6, color="rgb(0, 0, 0)"),
        hovertemplate="%{y:.0f} mg/dL, %{x|%H:%M:%p}",
    )

    # Delta BG Data
    delta_bg_trace = go.Scattergl(
        x=ji_df["cDateTime"],
        y=ji_df["deltaBg_smooth"],
        name="BG Velocity",
        showlegend=False,
        mode="markers",
        opacity=0.1,
        marker=dict(size=6, color="rgb(0, 0, 0)"),
        hovertemplate="%{y:.2f} mg/dL/min, %{x|%H:%M:%p}",
    )

    delta_bg_opt_lane_trace = go.Scattergl(
        x=ji_df["cDateTime"],
        y=ji_df["dG_optimal"],
        name="Optimal Settings",
        mode="markers",
        opacity=0.75,
        marker=dict(size=6, color="rgb(152, 78, 163)"),
        hovertemplate="%{y:.2f} mg/dL/min, %{x|%H:%M:%p}",
    )

    delta_bg_opt_1800_trace = go.Scattergl(
        x=ji_df["cDateTime"],
        y=ji_df["dG_nudge"],
        name="Nudged Settings",
        mode="markers",
        opacity=0.75,
        marker=dict(size=6, color="rgb(228, 26, 28)"),
        hovertemplate="%{y:.2f} mg/dL/min, %{x|%H:%M:%p}",
    )

    delta_bg_current_trace = go.Scattergl(
        x=ji_df["cDateTime"],
        y=ji_df["dG_current"],
        name="Current ISF & BR",
        mode="markers",
        opacity=0.75,
        marker=dict(size=6, color="rgb(255, 127, 0)"),
        hovertemplate="%{y:.2f} mg/dL/min, %{x|%H:%M:%p}",
    )

    cbg_trace_smooth = go.Scattergl(
        x=df["cDateTime"],  # cbg_ts["datetime_time"],
        y=df["mg_dL_smooth"],  # cbg_ts["value_mgdl"],
        name="All BGs Colored by<br>(Non-Insulin Velocity)",
        mode="markers",
        marker=dict(
            size=6,
            line=dict(width=0),
            color=df["dBg_non_insulin_" + insulin_model],
            colorscale="RdBu",
            cmin=-3,
            cmax=3,
        ),
        text=df["dBg_non_insulin_" + insulin_model],
        hovertemplate="%{y:.0f} mg/dL, %{text:.2f} mg/dL/min, %{x|%H:%M:%p}",
    )

    # BOLUS
    bolus_trace = go.Scattergl(
        name="Bolus Insulin",
        x=df["cDateTime"],
        y=df["normal"],
        mode="markers",
        marker=dict(symbol="triangle-down", color="lightskyblue", size=15,),
        hovertemplate="%{y:.2f}U, %{x|%H:%M:%p}",
    )

    # CARBS
    carb_trace = go.Scattergl(
        name="Carb Amount (g)",
        mode="markers + text",
        x=df.loc[df["grams"].notna(), "cDateTime"],
        y=df.loc[df["grams"].notna(), "normal"] + 5,
        marker=dict(color="gold", size=25),
        text=df.loc[df["grams"].notna(), "grams"].fillna(0),
        textposition="middle center",
        texttemplate="%{text:d}",
        hovertemplate="%{text:d}g, %{x|%H:%M:%p}",  # , g:%{text:d}g"
    )

    # Basal
    basal_trace = go.Scattergl(
        name="Delivered Basal Insulin (U/hr)",
        mode="lines",
        x=df["cDateTime"],
        y=df["basalRatePerMin"] * 60,
        line=dict(shape="vh", color="cornflowerblue"),
        fill="tozeroy",
        hovertemplate="%{y:.2f} U/hr, %{x|%H:%M:%p}",
    )

    basal_reg_trace = go.Scattergl(
        name="Scheduled Basal Rate (U/hr)",
        mode="lines",
        x=df["cDateTime"],
        y=df["scheduledBasalRate"],
        line=dict(shape="vh", color="black", dash="dash"),
        hovertemplate="%{y:.2f} U/hr, %{x|%H:%M:%p}",
    )

    # combine the plots
    bg_smooth_trace.yaxis = "y"
    delta_bg_trace.yaxis = "y2"
    delta_bg_opt_lane_trace.yaxis = "y2"
    delta_bg_opt_1800_trace.yaxis = "y2"
    delta_bg_current_trace.yaxis = "y2"
    cbg_trace_smooth.yaxis = "y3"
    bolus_trace.yaxis = "y4"
    carb_trace.yaxis = "y4"
    basal_trace.yaxis = "y5"
    basal_reg_trace.yaxis = "y5"

    data = [
        bg_smooth_trace,
        delta_bg_current_trace,
        delta_bg_opt_lane_trace,
        delta_bg_opt_1800_trace,
        delta_bg_trace,
        cbg_trace_smooth,
        carb_trace,
        bolus_trace,
        basal_trace,
        basal_reg_trace,
    ]

    # Set-up Layout
    layout = go.Layout(
        title=title,
        showlegend=True,
        yaxis=dict(title="BG<br>(mg/dL)", domain=[0.95, 1.0], range=[40, 400], fixedrange=True, tickvals=[40, 400],),
        yaxis2=dict(
            title="BG Velocity<br>(mg/dL/min)",
            domain=[0.65, 0.85],
            range=[-1, 1],
            tickvals=np.array([-snippet_threshold * 3, -snippet_threshold, snippet_threshold, snippet_threshold * 3]),
            fixedrange=True,
        ),
        yaxis3=dict(
            title="BG<br>(mg/dL)", domain=[0.4, 0.625], range=[40, 400], tickvals=[70, 180, 300], fixedrange=True,
        ),
        yaxis4=dict(title="Bolus<br>Events<br>(U)", domain=[0.175, 0.375], range=[0, 20], fixedrange=True,),
        yaxis5=dict(
            title="Basal<br>Insulin<br>(U/hr)",
            domain=[0, 0.15],
            range=[0, df["basalRatePerMin"].max() * 60 + 1],
            fixedrange=True,
        ),
        xaxis=dict(range=(START_DATE, END_DATE)),
        plot_bgcolor="#D3D3D3",
        dragmode="pan",
    )

    fig = go.Figure(data=data, layout=layout)

    return fig


def make_screen_results_fig(df):
    df.sort_values("rmse", inplace=True)
    df["opacity"] = df["rmse"].min() / df["rmse"]
    df["opacity"].fillna(np.min([0.15, df["opacity"].min() * 0.9]), inplace=True)
    df["rmse"].fillna("high", inplace=True)
    df.sort_index(inplace=True)

    fig = px.scatter(
        data_frame=df.round(decimals=2),
        x="isf",
        log_x=True,
        y="br",
        log_y=True,
        color="Legend",
        color_discrete_sequence=["rgba(51, 160, 44, 0.75)", "rgba(31, 120, 180, 0.75)", "rgba(255, 127, 0, 0.75)"],
        hover_data=["isf", "br", "rmse"],
        size="opacity",
        size_max=15,
    )

    title = "Preliminary search over ISF Range: 15 to 400 for userid {}<br>".format(
        userID
    ) + "Refining search to ISF Range: {} to {} mg/dL/U (shaded)".format(
        int(refined_isf_slice_min), int(refined_isf_slice_max)
    )

    layout = go.Layout(
        title=title,
        showlegend=True,
        xaxis=dict(title="ISF (mg/dL/U)", type="log", tickvals=isf_slice),
        yaxis=dict(title="Basal Rate (U/hr)", type="log", tickvals=br_slice[0 : len(br_slice) : 2],),
        plot_bgcolor="#D3D3D3",
        legend=dict(
            title="ISF & BR (Bigger = Better Fit)",
            orientation="v",
            xanchor="right",
            yanchor="top",
            x=0.99,
            y=0.985,
            bgcolor="rgba(255, 255, 255, 0.95)",
        ),
    )
    fig.update_layout(layout)

    fig.update_layout(
        shapes=[
            dict(
                name="refined search space",
                type="rect",
                x0=refined_isf_slice.min(),
                y0=refined_br_slice.max(),
                x1=refined_isf_slice.max(),
                y1=refined_br_slice.min(),
                fillcolor="rgba(25, 25, 25, 0.25)",
                layer="below",
                line_width=0,
            )
        ]
    )
    return fig


def unit_vector(v):
    return v / np.linalg.norm(v)


def get_nudged_results(current_isf, current_br, optimal_isf, optimal_br):
    v = np.array([optimal_isf - current_isf, optimal_br - current_br])
    unit_vector_v = unit_vector(v)

    # calculate the distance between the current settings and optimal settings
    normalized_distance = np.sqrt((((optimal_isf / current_isf) - 1) ** 2) + (((optimal_br / current_br) - 1) ** 2))

    ten_percent_normalized_distance = np.sqrt(2) * 0.10

    normalized_unit_vector_distance = np.sqrt(
        ((((current_isf + unit_vector_v[0]) / current_isf) - 1) ** 2)
        + ((((current_br + unit_vector_v[1]) / current_br) - 1) ** 2)
    )

    normalized_magnitude_10_percent = ten_percent_normalized_distance / normalized_unit_vector_distance

    if ten_percent_normalized_distance < normalized_distance:
        nudge_current_isf = int(current_isf + (unit_vector_v[0] * normalized_magnitude_10_percent))
        nudge_current_br = np.round(current_br + (unit_vector_v[1] * normalized_magnitude_10_percent), 2)

        # TODO: make this a test, nudged_normalized_distance == ten_percent_normalized_distance
        # nudged_normalized_distance = np.sqrt(
        #     (((nudge_current_isf / current_isf) - 1)**2)
        #     + (((nudge_current_br / current_br) - 1)**2)
        # )
    else:
        nudge_current_isf = optimal_isf
        nudge_current_br = optimal_br

    return nudge_current_isf, nudge_current_br


def make_refined_results_fig(df, title=""):
    nudge_current_isf = df.loc[df["Legend"].str.contains("Nudged"), "isf"].values[0]
    nudge_current_br = df.loc[df["Legend"].str.contains("Nudged"), "br"].values[0]
    df.sort_values("rmse", inplace=True)
    df["opacity"] = df["rmse"].min() / df["rmse"]
    df["opacity"].fillna(np.min([0.15, df["opacity"].min() * 0.9]), inplace=True)
    df["rmse"].fillna("high", inplace=True)
    df.sort_index(inplace=True)
    df = df.round(decimals=2)

    fig = px.scatter(
        data_frame=df,
        x="isf",
        log_x=True,
        y="br",
        log_y=True,
        color="Legend",
        color_discrete_sequence=[
            "rgba(51, 160, 44, 0.95)",
            "rgba(31, 120, 180, 0.95)",
            "rgba(255, 127, 0, 0.95)",
            "rgba(152, 78, 163, 1)",
            "rgba(227, 26, 28, 1)",
        ],
        hover_data=["isf", "br", "rmse"],
        size="opacity",
        size_max=15,
    )

    layout = go.Layout(
        title=title,
        showlegend=True,
        xaxis=dict(title="ISF (mg/dL/U)", type="log", tickvals=isf_slice),
        yaxis=dict(title="Basal Rate (U/hr)", type="log", tickvals=br_slice[0 : len(br_slice) : 2],),
        plot_bgcolor="#D3D3D3",
        legend=dict(
            title="ISF & BR with 50% Safety (shaded)",
            orientation="v",
            xanchor="right",
            yanchor="top",
            x=0.99,
            y=0.985,
            bgcolor="rgba(255, 255, 255, 0.95)",
        ),
    )

    fig.update_layout(layout)
    fig.update_layout(
        annotations=[
            dict(
                xref="x",
                yref="y",
                x=np.log10(nudge_current_isf),
                y=np.log10(nudge_current_br),
                axref="x",
                ayref="y",
                ax=np.log10(current_isf),
                ay=np.log10(current_br),
                showarrow=True,
                arrowcolor="rgba(255, 127, 0, 0.25)",
                arrowhead=4,
                arrowwidth=4,
            )
        ]
    )

    fig.update_layout(
        shapes=[
            dict(
                name="safety search space",
                type="rect",
                x0=safety_df["isf"].min(),
                y0=safety_df["br"].max(),
                x1=safety_df["isf"].max(),
                y1=safety_df["br"].min(),
                fillcolor="rgba(25, 25, 25, 0.25)",
                layer="below",
                line_width=0,
            )
        ]
    )
    return fig


# %% DOWNLOAD & PREPARE DATA
userid_of_shared_user = input("You acknowledge that this is exploratory (Press Return):\n")

if userid_of_shared_user in "":
    userid_of_shared_user = np.nan

weeks_of_data = np.int(input("How many weeks of data do you want to analyze? (2-4 is recommended)\n"))

current_isf = float(input("What is your current ISF while sleeping?\n"))
current_br = float(input("What is your current Basal Rate while sleeping?\n"))
ins_model_num = float(input("What Loop Insulin Model are you using?\n1=Adult, 2=Child, 3=Fiasp"))
if ins_model_num == 2:
    current_insulin_model = "adult"
elif ins_model_num == 3:
    current_insulin_model = "fiasp"
elif ins_model_num == 1:
    current_insulin_model = "adult"
else:
    print("we'll assume you meant Adult")
    current_insulin_model = "adult"
# TODO: add in some checks to make sure those values are correct

date_data_pulled = dt.datetime.now().strftime("%Y-%d-%mT%H-%M")

# data, responses = get_data_from_api(weeks_of_data=weeks_of_data, userid_of_shared_user=userid_of_shared_user)

# Get dataset from API
email = input("Enter the email address of your Tidepool account:\n")
if "bigdata" in email[:7]:
    password = os.environ.get("BIGDATA__PASSWORD")
else:
    password = getpass.getpass("Enter the password of your Tidepool account:\n")

print("\nGetting the last %d weeks of data..." % weeks_of_data)

data, dataset_userid = get_single_tidepool_dataset.get_dataset(
    weeks_of_data=weeks_of_data,
    userid_of_shared_user=userid_of_shared_user,
    email=email,
    password=password,
    return_raw_json=False,
)

print(len(data), "rows of data have been downloaded")

if pd.isna(userid_of_shared_user):
    userID = dataset_userid
else:
    userID = userid_of_shared_user

# %% CLEAN DATA
combined = clean_data(data)

# CALCULATE THE 3 LOOP INSULIN MODELS
combined = calc_loop_insulin_models(combined)

# RUN BRUTE FORCE METHOD OVER MESH GRID OF ISFs and BASAL RATES
combined = prep_insulin_snippets(combined)

# %% STEP 1: get a general sense of where solution is located using Lane criteria
snippet_threshold = 0.25
insulin_model = current_insulin_model

output_cols = ["Legend", "insulin_model", "rmse", "isf", "br"]
screen_lane_results = pd.DataFrame(columns=output_cols)
screen_1800_results = pd.DataFrame(columns=output_cols)

# lane's trivariate relationship (put link to spreadsheet here)
print("searching across ISF Range of 15 to 400...")
isf_slice = np.append(np.arange(15, 60, 5), np.append(np.arange(60, 120, 10), np.arange(120, 400, 50)),)
br_slice = np.round((((isf_slice / 494.45) ** (-1 / 0.7768)) / 24) / 0.025) * 0.025

for i, b in zip(isf_slice, br_slice):
    rule = "Lane's Rule"
    t_rmse = get_rmse_given_isf_br(i, b, combined, snippet_threshold, insulin_model)
    t_results = pd.DataFrame([[rule, insulin_model, t_rmse, i, b]], columns=output_cols)
    screen_lane_results = pd.concat([screen_lane_results, t_results], ignore_index=True)

br_slice_1800 = np.round((0.625 * 60 / isf_slice) / 0.025) * 0.025
for i, b in zip(isf_slice, br_slice_1800):
    rule = "1800 Rule"
    t_rmse = get_rmse_given_isf_br(i, b, combined, snippet_threshold, insulin_model)
    t_results = pd.DataFrame([[rule, insulin_model, t_rmse, i, b]], columns=output_cols)
    screen_1800_results = pd.concat([screen_1800_results, t_results], ignore_index=True)

# %% STEP 2: combine results and define a refined search space
rule = "Current (ISF={}, BR={})".format(int(current_isf), np.round(current_br, 2))
current_rmse = get_rmse_given_isf_br(current_isf, current_br, combined, snippet_threshold, insulin_model)
current_results = pd.DataFrame([[rule, insulin_model, current_rmse, current_isf, current_br]], columns=output_cols,)

screen_results = pd.concat([screen_lane_results, screen_1800_results, current_results], ignore_index=True,)
refined_isf_slice_min = int(screen_results.loc[screen_results["rmse"].notnull(), "isf"].min())
refined_isf_slice_max = int(screen_results.loc[screen_results["rmse"].notnull(), "isf"].max())
refined_isf_slice = np.arange(refined_isf_slice_min, refined_isf_slice_max + 1)
refined_br_slice = np.round((((refined_isf_slice / 494.45) ** (-1 / 0.7768)) / 24) / 0.025) * 0.025
refined_br_slice_1800 = np.round((0.625 * 60 / refined_isf_slice) / 0.025) * 0.025

# plot screening (preliminary) results
plot(make_screen_results_fig(screen_results))

# refined results
lane_results = pd.DataFrame(columns=output_cols)
rule1800_results = pd.DataFrame(columns=output_cols)

for i, b in zip(refined_isf_slice, refined_br_slice):
    rule = "Lane's Rule"
    t_rmse = get_rmse_given_isf_br(i, b, combined, snippet_threshold, insulin_model)
    t_results = pd.DataFrame([[rule, insulin_model, t_rmse, i, b]], columns=output_cols)
    lane_results = pd.concat([lane_results, t_results], ignore_index=True)

for i, b in zip(refined_isf_slice, refined_br_slice_1800):
    rule = "1800 Rule"
    t_rmse = get_rmse_given_isf_br(i, b, combined, snippet_threshold, insulin_model)
    t_results = pd.DataFrame([[rule, insulin_model, t_rmse, i, b]], columns=output_cols)
    rule1800_results = pd.concat([rule1800_results, t_results], ignore_index=True)

refined_results = pd.concat([lane_results, rule1800_results, current_results], ignore_index=True)

opt_indice = refined_results["rmse"].argmin()
optimal_isf = int(refined_results.loc[opt_indice, "isf"])
optimal_br = np.round(refined_results.loc[opt_indice, "br"], 2)
optimal_rmse = refined_results.loc[opt_indice, "rmse"]
rule = "Optimal (ISF={}, BR={})".format(optimal_isf, optimal_br)
opt_result = pd.DataFrame([[rule, insulin_model, optimal_rmse, optimal_isf, optimal_br]], columns=output_cols,)

just_insulin_df, avg_snippets_per_week, combined_temp = get_insulin_snippets_given_isf_br(
    combined, optimal_isf, optimal_br, snippet_threshold=snippet_threshold, insulin_model=insulin_model
)

# calculate safety parameters (i.e., br < optimal_br and isf > optimal_isf is safer)
optimal_egp = optimal_isf * optimal_br / 60
safety_df = pd.DataFrame(np.arange(optimal_isf, optimal_isf * 1.5, 1), columns=["isf"])
safety_df["br"] = optimal_egp * 60 / safety_df["isf"]
safety_df["Legend"] = "Safety Settings"
safety_df["insulin_model"] = "Adult"

for s_index in safety_df.index:
    safety_df.loc[s_index, "rmse"] = get_rmse_given_isf_br(
        safety_df.loc[s_index, "isf"], safety_df.loc[s_index, "br"], combined, snippet_threshold, insulin_model,
    )

nudge_current_isf, nudge_current_br = get_nudged_results(current_isf, current_br, optimal_isf, optimal_br)

rule = "Nudged 10% (ISF={}, BR={})".format(nudge_current_isf, nudge_current_br)
nudge_rmse = get_rmse_given_isf_br(nudge_current_isf, nudge_current_br, combined, snippet_threshold, insulin_model)
nudge_results = pd.DataFrame(
    [[rule, insulin_model, nudge_rmse, nudge_current_isf, nudge_current_br]], columns=output_cols
)


refined_results = pd.concat([refined_results, opt_result, nudge_results], ignore_index=True)
# nudge_current_isf = refined_results.loc[refined_results["Legend"].str.contains("Nudged"), "isf"].values[0]
# nudge_current_br = refined_results.loc[refined_results["Legend"].str.contains("Nudged"), "br"].values[0]

title = (
    "{} insulin snippets over {} weeks of {}'s data, <br>".format(
        int(avg_snippets_per_week * weeks_of_data), weeks_of_data, userID,
    )
    + "Optimal: ISF={}, BR={} ".format(optimal_isf, optimal_br,)
    + "Current Nudged by 10%: ISF={}, BR={}".format(nudge_current_isf, nudge_current_br)
)

plot(make_refined_results_fig(refined_results, title=title))


# %% make evidence plots
just_insulin_df["dG_current"] = get_delta_bg_from_isf_br_loop_insulin_curve(
    just_insulin_df["dG_insulin_loop_" + insulin_model].values, current_isf, current_br,
)

just_insulin_df["dG_optimal"] = get_delta_bg_from_isf_br_loop_insulin_curve(
    just_insulin_df["dG_insulin_loop_" + insulin_model].values, optimal_isf, optimal_br,
)

just_insulin_df["dG_nudge"] = get_delta_bg_from_isf_br_loop_insulin_curve(
    just_insulin_df["dG_insulin_loop_" + insulin_model].values, nudge_current_isf, nudge_current_br,
)

evidence_fig = make_evidence_plots(
    combined_temp, just_insulin_df, snippet_threshold, insulin_model, days_to_show=7, title=title,
)
plot(evidence_fig)
