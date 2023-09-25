import pandas as pd
from datetime import datetime, timedelta

def process_csv(input_file):
    # CSVファイルを読み込む
    df = pd.read_csv(input_file, header=None, dtype=str)  # dtype=strを指定して、すべてのデータを文字列として読み込む

    # 出力用のデータを保持するリストを初期化
    all_rows = []

    for i in range(0, len(df), 3):
        # 3行ずつデータを取得
        chunk = df.iloc[i:i+3]

        # 1列目は3行ずつの最初の行の値を使用
        first_column_value = chunk.iloc[0, 0]

        # avg_valuesリストを初期化し、1列目の値を追加
        avg_values = [first_column_value]

        # 2列目以降に対して平均を計算
        for col in range(1, len(chunk.columns)):
            # 列のデータを取得
            column_data = chunk.iloc[:, col]

            # 非数値データをNaNに変換
            column_data = pd.to_numeric(column_data, errors='coerce')

            # -9999を除外
            column_data = column_data[column_data != -9999]

            # 平均を計算: 全て-9999である場合は-9999を設定、それ以外は平均を計算
            avg_value = -9999 if column_data.isnull().all() else column_data.mean()
            avg_values.append(avg_value)

        # 出力用のデータリストに追加
        all_rows.append(avg_values)

    # 出力用のデータリストをDataFrameに変換
    output_df = pd.DataFrame(all_rows)

    # 平均値を同じCSVファイルに出力する（上書き）
    output_df.to_csv(input_file, index=False, header=False)  # index=Falseと

# 開始日と終了日を設定
start_date = datetime(1981, 1, 1)
end_date = datetime(2000, 12, 31)

# 日付の範囲でループ
current_date = start_date
while current_date <= end_date:
    input_filename = f'./../data/TYO/forcing/{current_date.strftime("%Y")}/TYO_{current_date.strftime("%Y")}_{current_date.strftime("%m_%d")}.csv'
    process_csv(input_filename)  # output_file パラメータは指定しない
    current_date += timedelta(days=1)
