/*
	Copyright 2023 Limeoats
	Copyright 2023 Samuel Ayala

   	Licensed under the Apache License, Version 2.0 (the "License");
   	you may not use this file except in compliance with the License.
   	You may obtain a copy of the License at

       	http://www.apache.org/licenses/LICENSE-2.0

   	Unless required by applicable law or agreed to in writing, software
   	distributed under the License is distributed on an "AS IS" BASIS,
   	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   	See the License for the specific language governing permissions and
   	limitations under the License.
*/

#pragma once

#include <string>
#include <optional>

namespace FlexGUI {
	const int max_path_length = 256;

	enum class FileDialogType {
		OpenFile,
		SaveFile,
		SelectFolder
	};
	enum class FileDialogSortOrder {
		Up,
		Down,
		None
	};

	class FileDialog {
		public:
			bool ready = false;
			bool file_dialog_open = false;
			FileDialogType type = FileDialogType::OpenFile;

			std::string file_dialog_current_path;

			void Show(char* buf);
			FileDialog();
	};
}
